

import numpy as np




f=open('gencode.v47.basic.annotation.gtf')
fw=open('ensembl_id_to_gene_symbols.csv','w')

fw.write('Ensembl_id,Gene_Symbol\n')

d={}

for line in f:
    l=line.split(';')

    if line.find('ENSG00000230427')!=-1:
        print(line)

    if len(l)>2:
        t=l[0].split('gene_id')

        answer=''
        for j in range(1,len(l)):
            if l[j].find('gene_name')!=-1:
                m=l[j].split('gene_name')
                answer=m[1]




        first = t[1].replace("'", "").replace('"', "")
        second = answer.replace("'", "").replace('"', "")

        name=first.strip()+','+second.strip()
        if name not in d:
            fw.write(name+'\n')
        d[name]=1
