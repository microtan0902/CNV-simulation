#-------------------------------------------------------------------------------
# Name:        simulation
# Purpose:	   Insert SNVs into fasta file
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

#demo command
# python simulation_insert_snp_to_fasta.py \
# -v /home/rjtan/data/NA12878_CEUtrio.HiSeq.WGS.b37_decoy.sort.chr1.phased.vcf \
# -r /home/rjtan/script/data/reference/hs37d5_chr1.fa \
# -c 1 \
# -s NA12878 \
# -o ./



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



# 67 server
# python simulation_insert_snp_to_fasta.py \
# -v /home/tanrenjie/data/CEUTrio/SNV/CEUtrio.WGS.b37_decoy_gatk3.3_HC_sort_PhaseByTransmission.vcf.gz \
# -r $ref \
# -c 1 \
# -s NA12878 \
# -o ./
parser = argparse.ArgumentParser(description='Aiming to insert SNPs into the reference genome.')
parser.add_argument('-v','--vcf',action="store", dest="vcf_file", required=True)
parser.add_argument('-r','--reference',action="store", dest="ref_file", required=True)
parser.add_argument('-c','--chromesome', action="store", dest="chromesome", required=False)
parser.add_argument('-s','--sample_name', action="store", dest="sample_name", required=True)
parser.add_argument('-o','--output',action="store", dest="output_path", required=True)

args = parser.parse_args()
vcf_file = args.vcf_file
ref_file = args.ref_file
chromesome_specific = args.chromesome
sample_name = args.sample_name
output_path = args.output_path
output_fasta_1 = (output_path+'/ref_'+sample_name+'_1.fa').replace('//','/')
output_fasta_2 = (output_path+'/ref_'+sample_name+'_2.fa').replace('//','/')
print "Output files:\t",output_fasta_1,output_fasta_2

#Read reference to memory
try:
	fasta_open = pysam.Fastafile(ref_file)
	print 'Fasta file open succeed.'
except:
	print 'Fasta file open failed.'
	exit(0)

try:
	length = fasta_open.get_reference_length(chromesome_specific)
	print length
	ref_revised_1_list = list(fasta_open.fetch(chromesome_specific,0,length))
	ref_revised_2_list = list(fasta_open.fetch(chromesome_specific,0,length))
	print 'Read %d bases of reference genome to memory.'%length
except:
	print 'Read reference genome to memory error.'

output_log_file1_name = output_fasta_1 + '.log'
output_log_file2_name = output_fasta_2 + '.log'
output_log_file1_fp = open(output_log_file1_name,'w')
output_log_file2_fp = open(output_log_file2_name,'w')
print output_log_file1_name,output_log_file2_name

output_log_file1_fp.write('chr\tpos\tref\talt\torigial\trevised\n')
output_log_file2_fp.write('chr\tpos\tref\talt\torigial\trevised\n')
#parse vcf
vcf_reader = vcf.Reader(open(vcf_file,'r'))
sample_name_list = vcf_reader.samples
print sample_name_list

if sample_name in sample_name_list:
	num = 0
	for record in vcf_reader:
		num += 1
		if record.is_snp:
			sample_GT = record.genotype(sample_name)['GT']
			# pdb.set_trace()
			if sample_GT != '0|0' and '0/0':
				#Revise both ref				 
				if sample_GT == '1|1' or sample_GT == '1/1':
					if record.REF == str(ref_revised_1_list[record.POS-1]).upper() and record.REF == str(ref_revised_2_list[record.POS-1]).upper():
						original_1 = str(ref_revised_1_list[record.POS-1]).upper()
						original_2 = str(ref_revised_2_list[record.POS-1]).upper()
						try:
							# print 'Revising both Refs %d\t%s:%d\t%s==%s,%s -> %s'%(num,record.CHROM,record.POS,record.REF,original_1,original_2,record.ALT)
							ref_revised_1_list[record.POS-1] = str(record.ALT[0]).upper()
							ref_revised_2_list[record.POS-1] = str(record.ALT[0]).upper()

							#Write to log file
							log1_str = record.CHROM+'\t'+str(record.POS)+'\t'+str(record.REF)+'\t'+str(record.ALT)+'\t'+str(original_1)+'\t'+str(ref_revised_1_list[record.POS-1])+'\n'
							log2_str = record.CHROM+'\t'+str(record.POS)+'\t'+str(record.REF)+'\t'+str(record.ALT)+'\t'+str(original_2)+'\t'+str(ref_revised_2_list[record.POS-1])+'\n'
							output_log_file1_fp.write(log1_str)
							output_log_file2_fp.write(log2_str)
						except:
							break
					else:
						print 'record.REF != reference. the coordinate may be wrong! please check.'
						break
				#Revise ref1
				elif sample_GT == '1|0':
					if record.REF == str(ref_revised_1_list[record.POS-1]).upper():
						original_1 = str(ref_revised_1_list[record.POS-1]).upper()
						try:
							# print 'Revising ref1 %d\t%s:%d\t%s==%s -> %s'%(num,record.CHROM,record.POS,record.REF,original_1,record.ALT)
							ref_revised_1_list[record.POS-1] = str(record.ALT[0]).upper()

							#Write to log file
							log1_str = record.CHROM+'\t'+str(record.POS)+'\t'+str(record.REF)+'\t'+str(record.ALT)+'\t'+str(original_1)+'\t'+str(ref_revised_1_list[record.POS-1])+'\n'
							output_log_file1_fp.write(log1_str)
						except:
							break
					else:
						print 'record.REF != reference. the coordinate may be wrong! please check.'
						break
				
				#Revise ref2
				elif sample_GT == '0|1':
					if record.REF == str(ref_revised_2_list[record.POS-1]).upper():
						original_2 = str(ref_revised_2_list[record.POS-1]).upper()
						try:
							# print 'Revising ref2 %d\t%s:%d\t%s==%s -> %s'%(num,record.CHROM,record.POS,record.REF,original_2,record.ALT)
							ref_revised_2_list[record.POS-1] = str(record.ALT[0]).upper()

							#Write to log file
							log2_str = record.CHROM+'\t'+str(record.POS)+'\t'+str(record.REF)+'\t'+str(record.ALT)+'\t'+str(original_2)+'\t'+str(ref_revised_2_list[record.POS-1])+'\n'
							output_log_file2_fp.write(log2_str)
						except:
							break
					else:
						print 'record.REF != reference. the coordinate may be wrong! please check.'
						break

				#in unphased condition. we revise a ref randomly
				else:
					ref_random = random.choice([1,2])
					if ref_random == 1:
						if record.REF == str(ref_revised_1_list[record.POS-1]).upper():
							original_1 = str(ref_revised_1_list[record.POS-1]).upper()
							try:
								# print 'Randomly, revising ref1 %d\t%s:%d\t%s==%s -> %s'%(num,record.CHROM,record.POS,record.REF,original_1,record.ALT)
								ref_revised_1_list[record.POS-1] = str(record.ALT[0]).upper()

								#Write to log file
								log1_str = record.CHROM+'\t'+str(record.POS)+'\t'+str(record.REF)+'\t'+str(record.ALT)+'\t'+str(original_1)+'\t'+str(ref_revised_1_list[record.POS-1])+'\n'
								output_log_file1_fp.write(log1_str)
							except:
								break
						else:
							print 'record.REF != reference. the coordinate may be wrong! please check.'
							break
					elif ref_random == 2:
						if record.REF == str(ref_revised_2_list[record.POS-1]).upper():
							original_2 = str(ref_revised_2_list[record.POS-1]).upper()
							try:
								# print 'Randomly, revising ref2 %d\t%s:%d\t%s==%s -> %s'%(num,record.CHROM,record.POS,record.REF,original_2,record.ALT)
								ref_revised_2_list[record.POS-1] = str(record.ALT[0]).upper()

								#Write to log file
								log2_str = record.CHROM+'\t'+str(record.POS)+'\t'+str(record.REF)+'\t'+str(record.ALT)+'\t'+str(original_2)+'\t'+str(ref_revised_2_list[record.POS-1])+'\n'
								output_log_file2_fp.write(log2_str)
							except:
								break
						else:
							print 'record.REF != reference. the coordinate may be wrong! please check.'
							break						

			else:
				print 'Genotype Error',sample_GT
				print record
				# exit(0)

	output_log_file1_fp.close()
	output_log_file2_fp.close()
	print('Finished to revise fasta file .')

	#Output the revised reference fasta
	print len(ref_revised_1_list),len(ref_revised_2_list)
	print output_fasta_1,output_fasta_2
	output_fasta(output_fasta_1,ref_revised_1_list,chromesome_specific)
	output_fasta(output_fasta_2,ref_revised_2_list,chromesome_specific)
else:
	print 'No %s in VCF file. Please check.'%sample_name
	