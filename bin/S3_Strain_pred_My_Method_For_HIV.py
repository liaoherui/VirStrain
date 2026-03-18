import re
import os
import getopt
import sys
import numpy as np
import json

opts,args=getopt.getopt(sys.argv[1:],"i:p:s:f:c:o:d:b:")
read_1=''
read_2=''
input_cls_file=''
out_file=''
id_map_file=''
input_kmr_file=''
db_dir=''
depth_cutoff=0.01

for opt,arg in opts:
	if opt=='-i':
		read_1=arg
	elif opt=='-p':
		read_2=arg
	elif opt=='-c':
		input_cls_file=arg
	elif opt=='-o':
		out_file=arg
	elif opt=='-d':
		id_map_file=arg
	elif opt=='-f':
		input_kmr_file=arg
	elif opt=='-b':
		db_dir=arg
dskc={} # Strain -> kmer number countdskc={} 
'''
Two dict used to visualize the depth figure
'''
dscl={} # Strain -> Column ->Label
dscf={} # Strain -> Column -> Freq
id2name={}
fi=open(id_map_file,'r')
while True:
	line=fi.readline().strip()
	if not line:break
	ele=line.split('\t')
	id2name[ele[0]]=ele[1]
	dscl[ele[0]]={}
	dscf[ele[0]]={}
	dskc[ele[0]]=0
id2cls={}
fc=open(input_cls_file,'r')
while True:
	line=fc.readline().strip()
	if not line:break
	if re.search('Cluster',line):
		cls=re.sub('>','',line)
	else:
		i=line.split('\t')[-1]
		id2cls[i]=cls
		
file_dir=sys.path[0]

#f=open(input_kmr_txt,'r')
dks={} #  kmer -> Strain  -> Column -> 0
dss={} # Strain -> map score
#dskc={} # Strain -> kmer number count
'''
Two dict used to visualize the depth figure
'''
#dscl={} # Strain -> Column ->Label
#dscf={} # Strain -> Column -> Freq
'''
while True:
	line=f.readline().strip()
	if not line:break
	ele=line.split('\t')
	dks[ele[0]]={}
	strain=re.split(',',ele[-1])
	column_raw=re.split(',',ele[1])
	column=[]
	for c in column_raw:
		rc=int(re.split('-',c)[0])
		column.append(rc)
	for s in strain:
		dks[ele[0]][s]={}
		dss[s]=0
		
		if s not in dscl:
			dscl[s]={}
		if s not in dscf:
			dscf[s]={}
		if s not in dskc:
			dskc[s]=0
	
		dskc[s]+=1
		#continue
		for c in column:
			#rc=int(re.split('-',c)[0])
			dscl[s][c]=len(strain)
			dscf[s][c]=0
			dks[ele[0]][s][c]=0
'''
dscl=json.load(open(db_dir+'/dscl.json'))
dscf=json.load(open(db_dir+'/dscf.json'))
dks=json.load(open(db_dir+'/dks.json'))
dss=json.load(open(db_dir+'/dss.json'))
dskc=json.load(open(db_dir+'/dskc.json'))
'''
json.dump(dscl,open("dscl.json",'w'))
json.dump(dscf,open("dscf.json",'w'))
json.dump(dks,open('dks.json','w'))
json.dump(dss,open('dss.json','w'))
json.dump(dskc,open('dskc.json','w'))
'''
print('Load Dict -> Done')
#exit()

if read_2=='':
	cmd1=file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if '+input_kmr_file+' -o Tem_VS.jf '+read_1
	cmd2=file_dir+'/jellyfish-linux dump -c Tem_VS.jf > Tem_Vs.fa'
	os.system(cmd1)
	os.system(cmd2)
else:
	cmd1=file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if '+input_kmr_file+' -o Tem_VS.jf '+read_1+' '+read_2
	cmd2=file_dir+'/jellyfish-linux dump -c Tem_VS.jf > Tem_Vs.fa'
	os.system(cmd1)
	os.system(cmd2)

fp=open('Tem_Vs.fa','r')
total=[]
while True:
	line=fp.readline().strip()
	if not line:break
	ele=line.split()
	if int(ele[-1])==0:continue
	total.append(int(ele[-1]))
min_depth=depth_cutoff*np.mean(total)

f2=open('Tem_Vs.fa','r')
dkc={}
dsvkc={} # Strain -> Valid kmer number count
#total=0
#ot=open('All_kmer_count.txt','w+')
while True:
	line=f2.readline().strip()
	if not line:break
	ele=line.split()
	if int(ele[-1])==0:continue
	#ss=','.join(dks[ele[0]])
	if int(ele[-1])<min_depth:continue
	#ot.write(line+'\t'+ss+'\n')
	dkc[ele[0]]=int(ele[-1])
	#total+=int(ele[-1])
	for s in dks[ele[0]]:
		if s not in dsvkc:
			dsvkc[s]={}
		#dsvkc[s]+=1
		
		#dss[s]+=float(ele[-1])/(float(len(dks[ele[0]])))
		dss[s]+=int(ele[-1])
		for c in dks[ele[0]][s]:
			dsvkc[s][c]=''
			dscf[s][c]+=int(ele[-1])
#print(len(dkc))
#exit()
dc={}
if os.path.exists(db_dir+'/head_file.txt'):
	fp=open(db_dir+'/head_file.txt','r')
	while True:
		line=fp.readline().strip()
		if not line:break
		ele=line.split('\t')
		dc[ele[0]]=ele[1]

o=open(out_file,'w+')
o.write('\t\tStrah_ID\tCls_info\tKmer_Hit_Num\tMap_Rate\tStrain_Depth\tStrain_info\n')
res=sorted(dss.items(),key=lambda d:d[1],reverse=True)
#o.write('Round-1\t'+res[0][0]+'\t'+str(res[0][1])+'\n')
depth_mp=float(res[0][1])/float(dskc[res[0][0]])
o.write('>>Most possible strains:\n')
if id2name[res[0][0]] not in dc:
	o.write('\t\t'+id2name[res[0][0]]+'\t'+id2cls[res[0][0]]+'\t'+str(res[0][1])+'\t'+str(len(dsvkc[res[0][0]]))+'/'+str(len(dscl[res[0][0]]))+'\t'+str(depth_mp)+'\tNA\n')
else:
	o.write('\t\t'+id2name[res[0][0]]+'\t'+id2cls[res[0][0]]+'\t'+str(res[0][1])+'\t'+str(len(dsvkc[res[0][0]]))+'/'+str(len(dscl[res[0][0]]))+'\t'+str(depth_mp)+'\t'+dc[id2name[res[0][0]]]+'\n')
mp_strain=id2name[res[0][0]]
mp_sid=res[0][0]
#exit()
cuid={res[0][0]:''}
'''
o3=open('check_used_kmer.txt','w+')
for k in dkc:
	if  res[0][0] in dks[k]:
		outs=','.join(dks[k])
		o3.write(k+'\t'+str(dkc[k])+'\t'+outs+'\n')
exit()
print(cuid)
'''
o.write('>>Other possible strains(Top5):\n')

c=0
#o2=open('check_rest_kmer.txt','w+')
op_strain=[]
op_sid=[]
while True:
	if c==5:break
	temc={}
	temc_knum={}
	#check=0
	count=0
	for k in dkc:
		check=0
		for cc in cuid:
			if cc in dks[k]:
				check=1
				#total=total-dkc[k]
				break
		if check==1:continue
		count+=1		
		#o2.write(k+'\t'+str(dkc[k])+'\t')
		#outs=','.join(dks[k])
		#o2.write(outs+'\n')
		for s in dks[k]:
			if s not in temc:
				temc[s]=dkc[k]
			else:
				temc[s]+=dkc[k]
			if s not in temc_knum:
				temc_knum[s]=1
				
			else:
				temc_knum[s]+=1
	#print(count)
	#exit()
	#if len(res)==0:break
	res=sorted(temc.items(),key=lambda d:d[1],reverse=True)
	if len(res)==0:break
	#exit()
	#print(res[:50])
	#exit()
	cuid[res[0][0]]=''
	op_sid.append(res[0][0])
	op_strain.append(id2name[res[0][0]])
	d=c+2
	cdepth=float(res[0][1])/float(temc_knum[res[0][0]])
	o.write('\t>>Left_Kmer_Num:'+str(count)+'\n')
	if id2name[res[0][0]] not in dc:
		o.write('\t\t'+str(id2name[res[0][0]])+'\t'+id2cls[res[0][0]]+'\t'+str(res[0][1])+'\t'+str(len(dsvkc[res[0][0]]))+'/'+str(len(dscl[res[0][0]]))+'\t'+str(cdepth)+'\tNA\n')
	else:
		o.write('\t\t'+str(id2name[res[0][0]])+'\t'+id2cls[res[0][0]]+'\t'+str(res[0][1])+'\t'+str(len(dsvkc[res[0][0]]))+'/'+str(len(dscl[res[0][0]]))+'\t'+str(cdepth)+'\t'+dc[id2name[res[0][0]]]+'\n')
	#o.write('Round-'+str(d)+'\t'+str(res[0][0])+'\t'+str(res[0][1])+'\n')
	c+=1
'''
res=sorted(dss.items(),key=lambda d:d[1],reverse=True)
for r in res[:10]:
	if not r[1]==0:
		o.write(r[0]+'\t'+str(r[1])+'\n')
'''
os.system('rm Tem_VS.jf Tem_Vs.fa')
print('Text report is done. Now will generate HTML report!')
ov1=open('Mps_ps_depth.csv','w+')
ov2=open('Ops_ps_depth.csv','w+')
ov1.write('ID,Column_ID,'+mp_strain+'_Freq,'+mp_strain+'_LNum\n')
carr=sorted(list(dscf[mp_sid].keys()))
i=1
for c in carr:
	ov1.write(str(i)+','+str(c)+','+str(dscf[mp_sid][c])+','+str(dscl[mp_sid][c])+'\n')
	i+=1

ov2.write('ID,Column_ID')
if len(op_strain)==0:
	ov2.write(',None,None\n')
else:
	#sid=0
	for s in op_strain:
		#sname=op_strain[sid]
		ov2.write(','+s+'_Freq'+','+s+'_LNum')
	ov2.write('\n')
	i=1
	for c in carr:
		ov2.write(str(i)+','+str(c))
		check=0
		check1=0
		check2=0
		for s in op_sid:
			if c not in dscf[s]:
				check+=1
				check1=1
			if c not in dscl[s]:
				check+=1
				check2=1
			if check==0:
				ov2.write(','+str(dscf[s][c])+','+str(dscl[s][c]))
			elif check==1:
				if check1==1:
					ov2.write(',0,'+str(dscl[s][c]))
				if check2==1:
					ov2.write(','+str(dscf[s][c])+',0')
			else:
				ov2.write(',0,0')
				
		ov2.write('\n')
		i+=1
		#sid+=1
