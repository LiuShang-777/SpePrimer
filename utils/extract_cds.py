# read fasta file and return fasta dictionary
def read_fa(cds_file):    
    with open(cds_file,'r') as file:
        name,seq,tmp=[],[],[]
        for line in file:
            line=line.strip()
            if line.startswith('>'):
                name.append(line[1:])
                seq.append(''.join(tmp))
                tmp=[]
            else:
                tmp.append(line)
        seq.append(''.join(tmp))
        seq=seq[1:]
        dic={}
        for i,j in zip(name,seq):
            dic[i]=j
    return dic




