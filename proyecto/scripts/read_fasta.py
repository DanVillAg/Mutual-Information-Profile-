import re
from Bio import SeqIO
from Mutual import generate_ami_profile, multiple_ami, random_sampling_corr, get_corr_against_gen
 
fasta1 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome1_human.fasta"
fasta2 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome7_human.fasta"
fasta3 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome14_human.fasta"
fasta4 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome17_human.fasta"

fasta5 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome1_mouse.fasta"
fasta6 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome7_mouse.fasta"
fasta7 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome14_mouse.fasta"
fasta8 = "/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome17_mouse.fasta"

fasta9='/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome_S_aureus.fasta'

fasta10='/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/chromosome1_c_elegans.fasta'
fasta11='/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/chromosome2_c_elegans.fasta'
fasta12='/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/chromosome3_c_elegans.fasta'
fasta13='/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/chromosome4_c_elegans.fasta'

fasta14='/home/dan/Escritorio/GenomicaComputacional/proyecto/data/raw_data/genome_E_coli.fasta'

# read FASTA file



def read_fasta(route_of_sequence, name_file=None):
    file = open(route_of_sequence, "r")
    for seq_record in SeqIO.parse(file, "fasta"):
        print(repr(seq_record.seq))
        genome = seq_record.seq.replace("N","")[0:500000]
        print(repr(genome))
        generate_ami_profile(genome, name=name_file)
        #print(len(genome))

def give_fasta_seq(sequence):
    '''
    Retorna unicamente la seccion de la secuencia para los genomas
    '''
    extracted = ''
    file = open(sequence, "r")
    for seq_record in SeqIO.parse(file, "fasta"):
        extracted = seq_record.seq.replace("N","")
    return extracted

if __name__=="__main__":

    h = 'human_'
    human_fasta= [fasta1,fasta2,fasta3,fasta4]
    #for idx, sequence in enumerate(human_fasta):
    #    read_fasta(sequence, name_file=h+str(idx))
    #filtered = [give_fasta_seq(gen)[0:500000] for gen in human_fasta]
    #multiple_ami(filtered, name=h+'all')

    m = 'mouse_'
    mouse_fasta= [fasta5,fasta6,fasta7,fasta8]

    #for idx, sequence in enumerate(mouse_fasta):
    #    read_fasta(sequence, name_file=m+str(idx))
    filtered = [give_fasta_seq(gen)[0:500000] for gen in mouse_fasta]
    #multiple_ami(filtered, name=m+'all')

    #get_corr_against_gen([give_fasta_seq(fasta9),give_fasta_seq(fasta14)])
    get_corr_against_gen([give_fasta_seq(fasta9)])


