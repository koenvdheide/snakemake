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
if (grepl('OS', fasta.header)) {
  pattern = "(?<=OS=).*(?=\\()"
} else {
  pattern = "(?<=\\[).*?(?=\\])"
}
species <- trimws(str_extract(string = fasta.header, pattern = pattern))
print(species)

######################################################################
# Determining paralogs and orthologs
######################################################################
species.splitted <- strsplit(species, ' ')[[1]]

# A function to check whether the name of the input species is different from 
# that of the entry in question. We cannot simply compare these as some
# Uniprot entries summarise strains (e.g. L.plantarum NC8/WCFS1).
is.paralog <- function(species.name) {
  pattern <- as.character(paste0("(?=.*", species.splitted,")", collapse=""))
  paralog <- grepl(pattern, species.name, perl = TRUE)
  return (paralog)
}

paralogs <- uniprot.table %>% 
  filter(is.paralog(.$Organism)) %>%
  droplevels()

orthologs <- uniprot.table %>%
  filter(!is.paralog(.$Organism)) %>%
  droplevels()

######################################################################
# Writing data
######################################################################
write.table(paralogs, row.names = F, quote = F, file = snakemake@output[[1]], sep = '\t')
write.table(orthologs, row.names = F, quote = F, file = snakemake@output[[2]], sep = '\t')













