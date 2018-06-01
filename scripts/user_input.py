import os.path
from Bio import SeqIO

def is_fasta(file_name):
    """
    This function tests whether the input file is in FASTA format using the
    BioPython SeqIO parser.
    :param file_name: Path to the file that should be checked
    :return: Whether the inputted file is in valid FASTA format
    """
    with open(file_name, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def is_valid(file_name):
    """
    This function checks whether the input file exists on the
    by the user specified location.
    :param file_name: The file path that should be checked
    :return: If the file exists or not
    """
    return os.path.isfile(file_name)

def copy_content(file_name):
    """
    This function copies the content of the by the user inputted file
    to a file accesible by snakemake. Hence, limiting errors due to
    file transfer or modifications to the file.
    :param file_name:
    :return:
    """
    with open(file_name) as in_file, open(snakemake.output[0],'w') as out_file:
        lines = in_file.readlines()
        out_file.writelines(lines)


'''
This script asks the user for a path to a FASTA file, whenever this file
does not exist or when it's not in FASTA format the user will be notified
and asked to give a valid file path.
'''
file_name = None
while (not file_name):
    file_name = input('Path to FASTA file: ')
    if file_name and is_valid(file_name) and is_fasta(file_name):
        copy_content(file_name)
    else:
        print("Specified file is not valid/FASTA")
        file_name = None
