#-------------------------------------------------------------------------------
# Name:        simulation
# Purpose:     Insert CNVs into a given reference genome.
#
# Author:      Renjie
#
# Created:     10/23/2017
# Copyright:   (c) Renjie 2017
# Licence:     MIT Licence
#-------------------------------------------------------------------------------
import vcf
import time
import pysam
import pdb
import argparse
import random
import os

def output_fasta(output_fasta_file,revised_fasta_list,chromesome_specific):
    try:
        output_fasta_fp = open(output_fasta_file,'w')
        revised_fasta_len = len(revised_fasta_list)
        output_fasta_fp.write('>'+chromesome_specific+'\n')
        for i in range(0,revised_fasta_len/60):
            output_fasta_fp.write(''.join(revised_fasta_list[i*60:i*60+60])+'\n')
        output_fasta_fp.write(''.join(revised_fasta_list[i*60+60:revised_fasta_len])+'\n')
        output_fasta_fp.close()
        print "The revised fasta file has been output to %s"%output_fasta_file
    except:
        print "output error"

def open_cnv_file(cnv_file):
    fp = open(cnv_file,'r')
    result_list = []
    for reader in fp:
        reader = reader.strip()
        one_line_list = reader.split('\t')
        
        result_list.append(one_line_list)
    return result_list

def insert_CNV(haplotype_str,hap_pointer_start,hap_pointer_stop):
    hap_pointer_stop = cnv_start - 1
    haplotype_str += fasta_open.fetch(cnv_chr,hap_pointer_start,hap_pointer_stop)
    print 'copy neutral region:',len(fasta_open.fetch(cnv_chr,hap_pointer_start,hap_pointer_stop))

    if cnv_GT == 'DEL':
        hap_pointer_start = cnv_stop + 1
        hap_pointer_stop = hap_pointer_start

    elif cnv_GT == 'DUP':
        hap_pointer_start = cnv_start
        hap_pointer_stop = cnv_stop
        dup_str = fasta_open.fetch(cnv_chr,hap_pointer_start,hap_pointer_stop)
        print 'DUP region:',len(dup_str)

        haplotype_str += dup_str + dup_str
        hap_pointer_start = cnv_start + cnv_length
        hap_pointer_stop = hap_pointer_start
        
    return [haplotype_str,hap_pointer_start,hap_pointer_stop]

# 67 server, demo command #############################
# python CNVsimulator.py \
# --fasta /home/tanrenjie/5_denovoCNV/simulation/ref_sim/ref_NA12878_1.fa \
# --cnv SIM_o_denovo.txt \
# -o ./
########################################################
parser = argparse.ArgumentParser(description='Aim to insert CNVs into a given reference genome.')
parser.add_argument('-f','--fasta',action="store", dest="fasta_file", required=True)
parser.add_argument('-chr','--chromesome', action="store", dest="chromesome", required=True)
parser.add_argument('-c','--cnv', action="store", dest="cnv_file", required=True)
parser.add_argument('-o','--output',action="store", dest="output_path", required=True)

args = parser.parse_args()
fasta_file = args.fasta_file
chromesome = args.chromesome
cnv_file = args.cnv_file
output_path = args.output_path

#prepare output path
output_hap1_file = output_path + '/SIM_temp_h1'+'/genome_rearranged.fasta'
output_hap2_file = output_path + '/SIM_temp_h2'+'/genome_rearranged.fasta'

if not os.path.exists(os.path.dirname(output_hap1_file)):
    os.makedirs(os.path.dirname(output_hap1_file))
    print 'path generated',os.path.dirname(output_hap1_file)

if not os.path.exists(os.path.dirname(output_hap2_file)):
    os.makedirs(os.path.dirname(output_hap2_file))
    print 'path generated',os.path.dirname(output_hap2_file)

#main function
try:
    fasta_open = pysam.Fastafile(fasta_file)
    print 'Fasta file open succeed.'
except:
    print 'Fasta file open failed.'
    exit(0)

length = fasta_open.get_reference_length(chromesome)
print "Total length:", length
print '--------------------------------------------------------------------'

num = 0
hap1_pointer_start = 0
hap1_pointer_stop = 0
hap2_pointer_start = 0
hap2_pointer_stop = 0

haplotype_1_str = ""
haplotype_2_str = ""

cnv_list = open_cnv_file(cnv_file)
for cnv_reader in cnv_list:
    num += 1
    cnv_chr = cnv_reader[0]
    cnv_start = int(cnv_reader[1])
    cnv_stop = int(cnv_reader[2])
    cnv_GT = cnv_reader[3]
    cnv_length = int(cnv_reader[4])
    cnv_allele_CN = cnv_reader[5]
    cnv_type = cnv_reader[6]

    print num, cnv_reader,hap1_pointer_start,hap1_pointer_stop,hap2_pointer_start,hap2_pointer_stop
    if cnv_allele_CN == '1|0':
        [haplotype_1_str,hap1_pointer_start,hap1_pointer_stop] = insert_CNV(haplotype_1_str,hap1_pointer_start,hap1_pointer_stop)

    elif cnv_allele_CN == '0|1':
        [haplotype_2_str,hap2_pointer_start,hap2_pointer_stop] = insert_CNV(haplotype_2_str,hap2_pointer_start,hap2_pointer_stop)

    elif cnv_allele_CN == '1|1':
        [haplotype_1_str,hap1_pointer_start,hap1_pointer_stop] = insert_CNV(haplotype_1_str,hap1_pointer_start,hap1_pointer_stop)
        [haplotype_2_str,hap2_pointer_start,hap2_pointer_stop] = insert_CNV(haplotype_2_str,hap2_pointer_start,hap2_pointer_stop)

    else:
        print "CNV allele error type: ", cnv_allele_CN

output_fasta(output_hap1_file,haplotype_1_str,chromesome)
output_fasta(output_hap2_file,haplotype_2_str,chromesome)

