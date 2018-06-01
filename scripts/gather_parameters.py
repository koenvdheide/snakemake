'''
As the user can change some parameters used in the pipeline it is important to include
these in the results, therefore these are saved whenever the pipeline is called.
'''

with open(snakemake.input[0]) as config, open(snakemake.output[0],'w') as out_file:
    header = []
    data = []
    for line in config:
        if ':' in line:
            # See config.yaml when wondering why the : split is required
            dat = line.strip().replace(' ','').split(':')
            header.append(dat[0])
            data.append(dat[1])
    out_file.write('\t'.join(header) + '\n')
    out_file.write('\t'.join(data) + '\n')
