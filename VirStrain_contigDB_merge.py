import re
import os,sys
import argparse

__author__="Liao Herui"

usage="VirStrain - An RNA virus strain-level identification tool for short reads."

def main():
    parser=argparse.ArgumentParser(prog='VirStrain_contigDB_merge.py',description=usage)
    parser.add_argument('-i','--merge_db',dest='merge_db',type=str,required=True,help="Dir of merged VirStrain database. eg: -i VirStrain_DB/SCOV2 or -i VirStrain_DB/SCOV2,VirStrain_DB/H1N1  --- Required")
    parser.add_argument('-o','--output_dir',dest='out_dir',type=str,help='Output dir of merged database. (default: current dir/VirStrain_DB_merge)')

    args=parser.parse_args()
    merge_db=args.merge_db
    out_dir=args.out_dir
    if not out_dir:
        out_dir='VirStrain_DB_merge'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    sp=re.split(',',merge_db)

    c=1
    d={}
    danno={}
    o=open(out_dir+'/merge_db.fa','w+')
    for s in sp:
        f=open(s+'/Pos-snp-kmer-all.fa','r')
        pre=re.split('/',s)[-1]
        while True:
            line=f.readline().strip()
            if not line:break
            if re.search('>',line):
                tem=line+'_'+pre
            else:
                if line not in d:
                    d[line]=1
                else:
                    d[line]+=1
                danno[line]=tem
        os.system('cp -rf '+s+' '+out_dir)
    for l in d:
        if d[l]>1:continue
        o.write(danno[l]+'\n'+l+'\n')


if __name__=='__main__':
    sys.exit(main())