suppressMessages(suppressWarnings(library(dplyr, quietly=TRUE)))
suppressMessages(suppressWarnings(library(stringr, quietly=TRUE)))

######################################################################
# Loading data
######################################################################
uniprot.table <- read.table(snakemake@input[[1]], sep = '\t', header = T, quote = "")
fasta.header <- readLines(snakemake@input[[2]])[1]

######################################################################
# Deriving species from input fasta
######################################################################
species <- trimws(str_extract(string = fasta.header, pattern = "(?<=OS=).*(?=\\()"))

######################################################################
# Determining paralogs and orthologs
######################################################################
paralogs <- uniprot.table %>% 
  filter(grepl(species,Organism)) %>%
  droplevels()

orthologs <- uniprot.table %>%
  filter(!grepl(species,Organism)) %>%
  droplevels()

######################################################################
# Writing data
######################################################################
write.table(paralogs, row.names = F, quote = F, file = snakemake@output[[1]], sep = '\t')
write.table(orthologs, row.names = F, quote = F, file = snakemake@output[[2]], sep = '\t')













