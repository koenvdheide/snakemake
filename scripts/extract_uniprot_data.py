import csv

NEWLINE = '\n'

def gene_names(entry, file):
    '''
    This function takes a gene entry and the file to which the name
    in the entry should be written to.
    :param entry: A UniProt gene entry
    :param file: A file to which the gene names should be written
    '''
    gene_names = entry.get('Gene.names',None)
    if gene_names:
        # only pick the first gene id as the rest is (occasionally) not mentioned in articles
        gene_name = gene_names.split()[0]
        file.write(gene_name + NEWLINE)

def gene_ids(entry, file):
    """
    This function takes a gene entry and the file to which the gene
    id in the entry should be written to.
    :param entry: A Uniprot entry
    :param file: A file to which the id should be written
    """
    gene_id = entry.get('Entry',None)
    if gene_id:
        file.write(gene_id + NEWLINE)

def kegg_ids(entry, file):
    '''
    This function takes a gene entry and a file to which the
    kegg ids should be written to.
    :param entry: A uniprot entry
    :param file: A file to which the keggs id should be written
    '''
    kegg_id = entry.get('Cross.reference..kegg.','').strip(';')
    if kegg_id:
        file.write(kegg_id + NEWLINE)


'''
This script reads Uniprot entries from an input file and 
subsquently extract relevant information form these
such as gene name, gene id and kegg id, which will 
all be written to separate files. 
'''
with    open(snakemake.input[0]) as in_file, \
        open(snakemake.output[0],'w') as gene_name_file, \
        open(snakemake.output[1],'w') as gene_id_file, \
        open(snakemake.output[2],'w') as kegg_id_file:
    reader = csv.DictReader(in_file, delimiter = '\t')
    for entry in reader:
        gene_names(entry, gene_name_file)
        gene_ids(entry, gene_id_file)
        kegg_ids(entry, kegg_id_file)
