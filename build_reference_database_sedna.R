###Reference library build/update, with taxonomy
###11/16/2022
###Code by Amy Van Cise and Megan Shaffer for MURI reference database management

### Load environment -----------------------------------------------------------

rm(list = ls())
library(rentrez)
library(tidyverse)
library(taxize)
library(Biostrings)


# navigate to an empty drive where you want the reference fasta files to go
# note that downstream analyses will read ALL fasta files in this location and use those to collect taxonomic information
reference_fasta_file_location <- ""
setwd(reference_fasta_file_location)

# create NCBI search terms
organism <- c("")
locus <- c("")

# create vector of terms to exclude
exclude <- c("shotgun", "predicted", "uncultured", "unclassified")

### Create reference fasta files -------------------------------------------------------------

for (i in 1:length(organism)) {
  for (j in 1:length(locus)){

search.term <- paste0("(", organism[i], "[ORGN] AND ", locus[j], ") ", "NOT (", paste(exclude, collapse = " OR "), ")")

ncbi.ids <- entrez_search(db="nucleotide", 
                      term=search.term,
                      use_history=TRUE)


fasta <- entrez_fetch(db="nucleotide", web_history=ncbi.ids$web_history, rettype="fasta")

locus <- locus %>% str_remove(pattern = "\\*")

write(fasta, file=paste(organism[i], "reference_fasta.fasta", sep = "_"))

  }
}

### Compile reference fastas ----------------------------------------------------

fasta.file.list <- list.files(., pattern = "*.fasta") 

fasta.files <- list()

for (i in 1:length(fasta.file.list)) {
fasta.files[[i]] <- readDNAStringSet(fasta.file.list[i])
fasta.files[[i]] <- fasta.files[[i]][!(grepl("UNVERIFIED", names(fasta.files[[i]])))]
}

for (i in 1:length(fasta.file.list)) {
  
### Generate species list for taxonomy grab ------------------------------------

fasta.og.names <- names(fasta.files[[i]]) %>% 
  as.data.frame(col.names = c("names")) %>% 
  mutate(.,ncbi_names = .) %>% 
  separate(.,col = .,into = c(NA,"Genus","species",NA,NA,NA), sep = " ") %>% 
  unite(full_species, Genus:species, sep = " ")

species.list <- fasta.og.names %>% 
  distinct(full_species) %>% 
  filter(!(grepl("mitochon", full_species))) %>% 
  filter(!(grepl("DNA", full_species))) %>% 
  filter(!(grepl("\\.", full_species))) %>% 
  pull(full_species)


### Get taxonomy information ---------------------------------------------------

tax.db <- classification(species.list, db = "itis", rows = 1) 
  
reference.db.taxonomy <- tax.db %>% 
  reshape::melt.list() %>% 
  select(-id, -value) %>% 
  pivot_wider(names_from = "rank", values_from = "name") %>% 
  select(L1, kingdom, phylum, class, order, family, genus, species) %>% 
  mutate(species = case_when(is.na(species) ~ L1, TRUE ~ species)) %>% 
  select(-genus) %>% 
  separate(species, into = c("genus", NA), remove = FALSE) %>% 
  relocate(genus, .before = species) %>% 
  select(-L1) %>% 
  mutate(kingdom = "Animalia") %>% 
  group_by(genus) %>% 
  fill(phylum, class, order, family, .direction = "downup") %>% 
  ungroup()

### Get common names for all species -------------------------------------------

common.names.list <- sci2comm(species.list, db = "itis", simplify = FALSE)
  
common.names <- do.call("rbind", common.names.list) %>% 
  rownames_to_column("species") %>% 
  separate(species, into = c("species", NA), sep = "\\.") %>% 
  filter(language == "English") %>% 
  group_by(tsn) %>% 
  slice_head() %>% 
  ungroup() %>% 
  select(-language, -tsn)

reference.db.taxonomy.common <- reference.db.taxonomy %>% 
  left_join(common.names, by = c("species" = "species"))

write.csv(reference.db.taxonomy.common, paste0(str_remove(fasta.file.list[i], ".fasta"), "_taxonomy.csv"))

### Match taxonomy to fasta files names ----------------------------------------

fasta.full.tax <- fasta.og.names %>% 
  left_join(reference.db.taxonomy.common, by = c("full_species" = "species")) %>% 
  mutate(species = full_species, .after = genus) %>% 
  unite(tax_name, kingdom:species, sep = ";") 

names(fasta.files[[i]]) <- fasta.full.tax[match(fasta.og.names$ncbi_names, fasta.full.tax$ncbi_names) , "tax_name"]  


### Write full taxonomy fasta file -----------------------------------------------
writeXStringSet(fasta.files[[i]], paste0(str_remove(fasta.file.list[i], ".fasta"), "_taxonomy.fasta"))
}
