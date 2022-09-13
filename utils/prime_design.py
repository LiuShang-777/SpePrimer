#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 15:45:22 2022

@author: ls
"""

import pandas as pd

def get_muscle_fa(input_muscle_fa):
    list_aln=[]
    with open(input_muscle_fa,'r') as file:
        dic_seq={}
        for line in file:
            line=line.strip()
            if ('MUSCLE' in line)|(len(line)==0):
                continue
            elif '*' in line:
                list_aln.append(line)
            else:
                line=line.strip().split(' ')
                if line[0] in dic_seq.keys():
                    dic_seq[line[0]].append(line[-1])
                else:
                    dic_seq[line[0]]=[line[-1]]
    return dic_seq

def trans_seq_to_dat(dic_seq):
    table=str.maketrans('ATCG-','12340') 
    for i in dic_seq.keys():
        dic_seq[i]=''.join(dic_seq[i]).translate(table)
        dic_seq[i]=[int(z) for z in dic_seq[i]]
    dat_seq=pd.DataFrame(dic_seq)
    return dat_seq

def judge_specific(dat_result,index_,treat_index):    
    test=dat_result.loc[dat_result.index.isin([index_,treat_index])]
    test1=test.sum(axis=0).min(axis=0)
    return test1

def trans_num_to_str(s1):
    table=str.maketrans('12340','ATCG-')
    s1=s1.translate(table)
    return s1
    
def reverse_str(s1):
    table=str.maketrans('ATCG-','TAGC-') 
    s1=s1.translate(table)
    s1=s1[::-1]
    return s1

def get_seq(start,end,list_seq,length):
    start_s=[]
    count=0
    for i in list_seq[start:]:
        if count<length:    
            if i!=0:
                count+=1
                start_s.append(str(i))
    count=0
    end_s=[]
    for i in list_seq[end-length+1:]:
        if count<length:
            if i!=0:
                count+=1
                end_s.append(str(i))
    return (''.join(start_s),''.join(end_s))
            
def tm_cal(seq):
    if len(seq)<=25:
        return 2*(seq.count('A')+seq.count('T'))+4*(seq.count('C')+seq.count('G'))
    else:
        return 81.5+0.41*((seq.count('C')+seq.count('G'))/len(seq))-600/len(seq)

def judge_tandem(seq):
    record=[0]
    for i,j in zip(seq[:-1],seq[1:]):
        if i==j:
            record.append(record[-1]+1)
        else:
            record.append(0)
    record.sort()
    return '%.2f'%(100*(record[-1]+1)/len(seq))

def judge_complementary(seq1,seq2):
    table=str.maketrans('ATCG-','TAGC-')
    seq2_=seq2.translate(table)
    cnt=0
    for i,j in zip(seq1,seq2_):
        if i==j:
            cnt+=1
        else:
            continue
    return '%.2f'%(100*cnt/len(seq1))

def implement(index_start,index_end,seq_list,length):    
    test_start,test_end=get_seq(index_start,index_end,seq_list,length)
    test_start=trans_num_to_str(test_start)
    test_end=trans_num_to_str(test_end)
    test_end=reverse_str(test_end)       
    tm_start=tm_cal(test_start)
    tm_end=tm_cal(test_end)
    tandem_rate_start=judge_tandem(test_start)
    tandem_rate_end=judge_tandem(test_end)
    inter_complementary=judge_complementary(test_start,test_end)
    start_complementary=judge_complementary(test_start,test_start[::-1])
    end_complementary=judge_complementary(test_end,test_end[::-1])
    return ([test_start,test_end,tm_start,tm_end,tandem_rate_start,tandem_rate_end,
             inter_complementary,start_complementary,end_complementary])        