suppressMessages(suppressWarnings(library(dplyr, quietly=TRUE)))
suppressMessages(suppressWarnings(library(stringr, quietly=TRUE)))

######################################################################
# Loading data
######################################################################
# - reading the header specifications (and cleaning it)
cat('[INFO] extracting header from BLAST output\n')
blast.results <- snakemake@input[[1]]
header <- scan(blast.results, skip = 3, nlines = 1, what = list('character'), sep = ',')
header <- gsub('# Fields:|\\s','',header[[1]])

# - reading data and attaching header
cat('[INFO] Reading matches\n')
data.table <- read.table(blast.results, sep = '\t', header = T)
colnames(data.table) <- header

######################################################################
# Filtering matches
######################################################################
cut.off <- snakemake@params[[1]]
id.cut.off <- snakemake@params[[2]]
matches <- data.table %>%
  filter(evalue < cut.off) %>%
  filter(`%identity` > id.cut.off)
cat(paste0('[INFO] Found ', nrow(matches), ' matches fitting the filtering criteria\n'))

######################################################################
# Saving significant matches
######################################################################
cat('[INFO] Saving significant matches\n')
write.table(matches, file = snakemake@output[[1]], sep = '\t',
            row.names = F)

######################################################################
# Saving significant ids
######################################################################
cat('[INFO] Saving significant IDs\n')
match.ids <- unlist(lapply(as.character(matches$subjectacc.ver), FUN = function(x) {unlist(strsplit(x,":"))[2]}))
write.table(data.frame(match.ids), file = snakemake@output[[2]], col.names = F,
            row.names = F, quote = F)
