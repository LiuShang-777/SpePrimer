#get input blast list; output name list
def get_name(input_blast_list):
    dic,per={},[]
    for line in input_blast_list:
        line=line.strip().split('\t')
        dic[line[1]]=float(line[2])
        if int(float(line[2]))!=100:
            per.append(float(line[2]))
    per.sort()
    result=[]
    for i in dic.keys():
        if dic[i] in per:
            result.append(i)
    return result







