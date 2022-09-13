#get blast file and output target blast list
def get_blast(input_blast,name):
    with open(input_blast,'r') as file:
        list_result=[]
        for line in file:
            line=line.strip().split('\t')
            if line[0]==name:
                list_result.append('\t'.join(line))
    return list_result




