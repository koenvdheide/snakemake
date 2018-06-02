"""
Author: Rick Beeloo
Aim: A simple worfklow to aid in unraveling the function of a protein
Date: 09-05-2018
Run: snakemake -s Snakefile
Latest modification:
    - removed extract scripts and used command line commands as replacement
"""

configfile: "config.yaml"

rule user_input:
    message:
        '--- File input ---'
    output:
        'scripts/data/sequence.fasta'
    script:
        'scripts/user_input.py'

rule save_parameters:
    message:
        '--- Saving the used parameters ---'
    input:
        'config.yaml'
    output:
        'scripts/data/parameters.txt'
    script:
        'scripts/gather_parameters.py'

rule phobius:
    message:
        "--- Performing Location prediction (Phobius) ---"
    input:
        'scripts/data/sequence.fasta'
    output:
        'scripts/data/phobius_result.visual-png.png'
    shell:
        'perl scripts/phobius_lwp.pl --email {config[mail]}' +
        ' --outfile scripts/data/phobius_result {input} --format grp'

rule interpro:
    message:
        "--- Performing InterProScan ---"
    input:
        'scripts/data/sequence.fasta'
    output:
        'scripts/data/proscan_result.svg.svg'
    shell:
        'scripts/iprscan5_urllib3.py --email {config[mail]}' +
        ' --outfile scripts/data/proscan_result --sequence {input}'

rule blast:
    message:
        "--- Performing BLAST ---"
    input:
        'scripts/data/sequence.fasta'
    params:
        blast_database = config['blast_database'],
        max_align = config['blast_max_alignments']
    output:
        'scripts/data/blast_results.out.txt',
        'scripts/data/blast_results.complete-visual-svg.svg'
    shell:
        'perl scripts/ncbiblast_lwp.pl --email {config[mail]} --sequence {input}' +
        ' --database {params.blast_database} --stype protein --program blastp --align 7' +
        ' --alignments {params.max_align} --outfile scripts/data/blast_results'

rule analyze_matches:
    message:
        "--- Filtering BLAST matches ---"
    input:
        'scripts/data/blast_results.out.txt'
    params:
        cut_off = config['e_value_cut_off'],
        identity = config['identity_cut_off']
    output:
        'scripts/data/significant_matches.txt',
        'scripts/data/significant_gene_ids.txt'
    script:
        'scripts/filter_matches.R'

rule uniprot_annotation:
    message:
        "--- Gathering uniprot records ---"
    input:
        'scripts/data/significant_gene_ids.txt'
    output:
        'scripts/data/uniprot_records.txt'
    shell:
        "query=$(python3 scripts/formatter.py {input} '+OR+' q)\n"
        'wget "http://www.uniprot.org/uniprot/?query=$query&format=tab&columns=id,entry%20name,protein%20names,ec,organism,organism-id,genes,comment(PATHWAY),go(biological process),database(kegg)" --output-document {output}'

rule orthologs:
    message:
        '--- Splitting orthologs and paralogs ---'
    input:
        'scripts/data/uniprot_records.txt',
        'scripts/data/sequence.fasta'
    output:
        'scripts/data/uniprot_para_records.txt',
        'scripts/data/uniprot_ortho_records.txt'
    script:
        'scripts/paralogs_orthologs.R'

rule uniprot_sequences:
    message:
        "--- Gathering protein sequences ---"
    input:
        'scripts/data/uniprot_ids.txt'
    output:
        'scripts/data/uniprot_sequences.fasta'
    shell:
        "query=$(python3 scripts/formatter.py {input} '+OR+' q)\n"
        'wget "http://www.uniprot.org/uniprot/?query=$query&format=fasta" --output-document {output}'

rule extract_data:
    message:
        "--- Extracting sequences and gene names from UnipProt records ---"
    input:
        'scripts/data/uniprot_ortho_records.txt'
    output:
        'scripts/data/uniprot_genes.txt',
        'scripts/data/uniprot_ids.txt',
        'scripts/data/kegg_ids.txt'
    script:
        'scripts/extract_uniprot_data.py'

rule msa:
    message:
        "--- Performing MSA ---"
    input:
        'scripts/data/uniprot_sequences.fasta'
    output:
        'scripts/data/t_coffee_msa.aln-clustalw.clustalw',
        'scripts/data/t_coffee_msa.phylotree.ph'
    shell:
        'scripts/tcoffee_lwp.pl --email {config[mail]} --outfile scripts/data/t_coffee_msa {input}'

rule clustal_to_fasta:
    message:
        '--- Converting clustalW to FASTA format ---'
    input:
        'scripts/data/t_coffee_msa.aln-clustalw.clustalw'
    output:
        'scripts/data/fasta_msa.out.txt'
    shell:
        'scripts/emboss_seqret_lwp.pl {input} --inputformat clustal' +
        ' --outputformat fasta --email {config[mail]} --stype protein --outfile scripts/data/fasta_msa'

rule literature:
    message:
        '--- Querying PubMed --- '
    input:
        'scripts/data/uniprot_genes.txt'
    params:
        max_articles = config['maximum_articles_retrieve']
    output:
        'scripts/data/literature_results.txt'
    shell:
        "query=$(python3 scripts/formatter.py {input} ' OR ')\n"
        'esearch -db gene -query \"$query\" | elink -target pubmed | efetch -format xml --stop {params.max_articles}' +
        ' | xtract -pattern PubmedArticle -element MedlineCitation/PMID ArticleTitle PubDate/Year Journal/ISSN  Abstract/AbstractText> {output}'

rule interactions:
    message:
        '--- Getting interactions from STRING ---'
    input:
        'scripts/data/uniprot_genes.txt'
    params:
        max_con = config['max_numb_con']
    output:
        'scripts/data/network.png'
    shell:
        "query=$(python3 scripts/formatter.py {input} '%0d')\n"
        'wget "https://string-db.org/api/image/network?identifiers=$query&limit={params.max_con}" --output-document {output}'

rule kegg:
    message:
        '--- Getting pathways from KEGG ---'
    input:
        'scripts/data/kegg_ids.txt'
    output:
        'scripts/data/kegg_table.txt'
    shell:
        "query=$(python3 scripts/formatter.py {input} ',')\n"
        "curl http://togows.org/entry/kegg-genes/$query.json" +
        "| cat | jq -r '.[] | .genes_id as $id | .pathways | to_entries[] | " +
        "[$id, .key, .value] | @tsv' > {output}"


rule results:
    message:
        '--- Generating result page  ---'
    input:
        'scripts/data/phobius_result.visual-png.png',
        'scripts/data/parameters.txt',
        'scripts/data/proscan_result.svg.svg',
        'scripts/data/blast_results.complete-visual-svg.svg',
        'scripts/data/uniprot_records.txt',
        'scripts/data/t_coffee_msa.phylotree.ph',
        'scripts/data/fasta_msa.out.txt',
        'scripts/data/literature_results.txt',
        'scripts/data/network.png',
        'scripts/data/kegg_table.txt'
    shell:
        "R -e \"library(rmarkdown); render('scripts/results.Rmd', output_dir = 'scripts/data')\""
