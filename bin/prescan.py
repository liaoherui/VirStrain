import re
import os

def scan(ingenome,db):
    file_dir=os.path.split(os.path.abspath(__file__))[0]
    fd=open(db,'r')
    d={}
    while True:
        line=fd.readline().strip()
        if not line:break
        if re.search('>',line):
            pre=re.split('_',line)
            pre='_'.join(pre[1:])
        else:
            d[line]=pre
    #print(file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if '+db+' -o Tem_Vs.jf '+ingenome)
    os.system(file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if '+db+' -o Tem_Vs.jf '+ingenome)
    #print(file_dir+'/jellyfish-linux dump  Tem_Vs.jf > Tem_Vs.fa')
    os.system(file_dir+'/jellyfish-linux dump -c Tem_Vs.jf > Tem_Vs.fa') 
    f=open('Tem_Vs.fa','r')
    dr={}
    while True:
        line=f.readline().strip()
        if not line:break
        ele=line.split()
        if len(ele)<2:continue
        if int(ele[-1])>0:
            if d[ele[0]] not in dr:
                dr[d[ele[0]]]=0
            dr[d[ele[0]]]+=1
    #print(dr)
    res=sorted(dr.items(),key=lambda d:d[1],reverse=True)
    if len(res)<1:
        return 'NA',0
    tsp=res[0][0]
    #print(tsp)
    #exit()
    return tsp,res[0][1]
