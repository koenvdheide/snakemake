
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('scripts/data/uniprot_records.txt', 'scripts/data/sequence.fasta'),
    output = list('scripts/data/uniprot_para_records.txt', 'scripts/data/uniprot_ortho_records.txt'),
    params = list(),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("identity_cut_off" = 60, "max_numb_con" = 20, "e_value_cut_off" = 0.05, "blast_database" = 'uniprotkb_bacteria', "min_con_score" = 400, "blast_max_alignments" = 250, "mail" = 'r.beeloo@outlook.com', "maximum_articles_retrieve" = 10),
    rule = 'orthologs'
)
######## Original script #########
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













