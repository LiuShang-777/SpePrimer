#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:31:01 2022

@author: ls
"""
import sys
import os
import pandas as pd
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-gene",help="gene name of target")
parser.add_argument("-cds",help="CDS of the whole genome")
parser.add_argument("-blast",help="the whole genome blast result")
parser.add_argument("-out",help="the output results of primer design")
parser.add_argument("-maxlen",help="maximum length of primer products")
parser.add_argument("-minlen",help="minimum length of primer products")
parser.add_argument("-primerlen",help="primer length")
args=parser.parse_args()

sf_path=sys.path[0]
sys.path.append(sf_path)
from utils import extract_cds, get_names, extract_blast, prime_design

blast_file=args.blast
cds_file=args.cds
gene=args.gene
distance_min=int(args.minlen)
distance_max=int(args.maxlen)
prime_length=int(args.primerlen)
output_file=args.out

'''
# args parser
blast_file='/home/ls/primer/database/zm12.blast'
cds_file='/home/ls/primer/database/zm12.cds.fasta'
gene='CRI12_D13G1278.1'
distance_min=200
distance_max=300
prime_length=20
output_file='/home/ls/primer/database/CRI12_D13G1278.output'
'''

print("Primer design for %s is starting\n"%gene)
#reading blast file
print("Reading the blast file now\n")
all_blast_records=extract_blast.get_blast(blast_file,gene)
#homo_genes=get_names.get_name(all_blast_records)
homo_genes=[i.split('\t')[1] for i in all_blast_records]
homo_genes=list(set(homo_genes))
#reading cds file
print("Reading the cds file now\n")
all_cds=extract_cds.read_fa(cds_file)
gene_seq_dic={}
for gene_ in homo_genes:
    gene_seq_dic[gene_]=all_cds[gene_]
if os.path.exists(sf_path+'/tmp/'):
    pass
else:    
    os.makedirs(sf_path+'/tmp/')
with open(sf_path+'/tmp/'+'%s.cds.fasta'%gene,'w') as file:
    for line in gene_seq_dic.keys():
        file.write('>'+line+'\n')
        file.write(gene_seq_dic[line]+'\n')
    
#performing multiple alignments 
print("Performing multi-alignments on %d sequences\n"%len(homo_genes))
os.system("muscle -in %s/tmp/%s.cds.fasta -out %s/tmp/%s.aln.fasta -clw"%(sf_path,gene,sf_path,gene))

#start designing specifical primers
print("Starting specific primer design for %s\n"%gene)
aln_fasta_file=sf_path+'/tmp/'+'%s.aln.fasta'%gene
dic_target_fasta=prime_design.get_muscle_fa(aln_fasta_file)
#using numeric sequences to detect unalign sites
dat_target_fasta=prime_design.trans_seq_to_dat(dic_target_fasta)
dat_no_aln=dat_target_fasta.loc[dat_target_fasta[homo_genes].mean(axis=1)!=dat_target_fasta[homo_genes[0]]]
dat_no_aln_result=pd.DataFrame()
for i in homo_genes:
    if i !=gene:    
        dat_no_aln_result[i]=abs(dat_no_aln[i]-dat_no_aln[gene])

#select loc pairs of proper length
print("Judging the specifity of potential primers")
index_list=dat_no_aln_result.index.tolist()
record_list=[[0,0]]
cnt=0
for i in index_list:
    treat_index=[index_ for index_ in index_list if index_>=i+distance_min]
    for index_ in treat_index:
        if prime_design.judge_specific(dat_no_aln_result,i,index_)==0:
            continue
        else:
            if (i+1) >=20+record_list[-1][0]:    
                record_list.append([i+1,index_+1])
record_list=record_list[1:]

#record primers
test={}
dic_primer_product={}
for index_pair in record_list:
    index_range=range(index_pair[0],index_pair[1])
    dat_product=dat_target_fasta.loc[dat_target_fasta.index.isin(index_range)]
    dat_product=dat_product.sort_index()
    tmp=prime_design.implement(index_pair[0],index_pair[1],dat_target_fasta[gene],prime_length)
    if (len(tmp[0])==prime_length)&(len(tmp[1])==prime_length):
        idx=str(index_pair[0])+'_'+str(index_pair[1])
        test[idx]='\t'.join([str(i) for i in tmp])
        align_=dat_product.loc[dat_product.mean(axis=1)==dat_product[homo_genes[0]]]
        dic_primer_product[idx]=align_.shape[0]/dat_product.shape[0]

#relocation of primer in target sequences
target_sequence=prime_design.trans_num_to_str(''.join([str(i) for i in dat_target_fasta[gene]]))
target_sequence=''.join([i for i in target_sequence if i!='-'])
os.system("rm -rf %s/tmp"%sf_path)
#geting information recorded in the output file
products=[]
for i in test.keys():
    test[i]=test[i].split('\t')
    start_site=target_sequence.find(test[i][0])
    end_site=target_sequence.find(prime_design.reverse_str(test[i][1]))
    CG_F=str(round((test[i][0].count('C')+test[i][0].count('G'))/len(test[i][0]),2))
    CG_R=str(round((test[i][1].count('C')+test[i][1].count('G'))/len(test[i][1]),2))
    if ((end_site-start_site)<=distance_max)&((end_site-start_site)>=distance_min):
        products.append('\t'.join(test[i])+'\t'+str(round(dic_primer_product[i],2))+'\t'+CG_F+'\t'+CG_R+'\t'+str(start_site)+'\t'+str(end_site)+'\t'+str(end_site-start_site))
with open(output_file,'w') as file:
    file.write('\t'.join(['forward','reverse','tm_fr','tm_rv','tandem_rate_fr','tandem_rate_rv','inter_com','self_com_fr','self_com_rv','pd_similarity','CG_F','CG_R','start','end','length'])+'\n')
    for i in products:
        file.write(i+'\n')    

