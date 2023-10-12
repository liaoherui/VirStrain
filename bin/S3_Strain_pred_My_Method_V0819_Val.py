import re
import os
import getopt
import sys
import subprocess
#from collections import Counter
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
import numpy as np
#import pandas as pd
### Dehan Method Check
#from scipy import sparse
#from sklearn.model_selection import ShuffleSplit
#from sklearn.linear_model import Lasso, LassoCV
#from scipy.stats import pearsonr

opts,args=getopt.getopt(sys.argv[1:],"i:p:k:m:s:f:c:b:d:r:o:")
read_1=''	# SE or long reads
read_2=''	# PE reads - read set 2
snp_kmr_file='' #kmer -> snp_pos file
snp_kmr_fa='' # kmr fasta
matrix_file='' # Sno-pos matrix file
cls_file=''
sub_kmr_file=''
out_dir=''
rks=''
db_dir=''

k=25		# Default k size=25
#min_depth_rate=0.1
min_depth_percentile=10
max_depth_percentile=90
min_depth_absolute=2
#min_depth_rate=0.05
min_depth_rate=0.05
#threads=4
for opt,arg in opts:
	if opt=='-i':
		read_1=arg
	elif opt=='-p':
		read_2=arg
	elif opt=='-k':
		k=int(arg)
	elif opt=='-s':
		snp_kmr_file=arg
	elif opt=='-m':
		matrix_file=arg
	elif opt=='-f':
		snp_kmr_fa=arg
	elif opt=='-c':
		cls_file=arg
	elif opt=='-b':
		sub_kmr_file=arg
	elif opt=='-d':
		db_dir=arg
	elif opt=='-r':
		rks=int(arg)
	elif opt=='-o':
		out_dir=arg


## Static variables
BASE_ORDER=['A','T','G','C']
BASE_P = {'A': [1, 0, 0, 0],'C':[0,1,0,0],'G':[0,0,1,0],'T':[0,0,0,1],}
'''
CV_NITER = 20
NALPHA = 50
MAX_NITER = 5000
TEST_SIZE = 0.5
'''
file_dir=sys.path[0]
#print(file_dir)
#exit()


######### Step-1 Load Pre-build File to memory ####
## Kmer -> POS-SNP
f1=open(snp_kmr_file,'r')
dkps={}  # kmr -> {pos-snp:1,......}
pos=[]
dpsc={} # pos-snp:  num
while True:
	line=f1.readline().strip()
	if not line:break
	ele=line.split('\t')
	dkps[ele[0]]=''
	ps=re.split(',',ele[1])
	for e in ps:
		dkps[ele[0]]=e # Set to 1 for Counter -> Dict Merge
		dpsc[e]=0

## Build pos-snp freq array

f3=open(matrix_file,'r')
fl=f3.readline().strip()
pos_snp=re.split('\t',fl) # Head line arr
#print(np.where(np.array(pos_snp)=='8946-T')[0][0])

# Run jellyfish to get kmer counting result
if read_2=='':
	if re.split('\.',read_1)[-1]=='gz':
		cmd1='zcat '+read_1+' | '+file_dir+'/jellyfish-linux count /dev/fd/0 -m 25 -s 100M -t 8 --if '+snp_kmr_fa+' -o Tem_VS.jf '
	else:
		cmd1=file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if '+snp_kmr_fa+' -o Tem_VS.jf '+read_1
	cmd2=file_dir+'/jellyfish-linux dump -c Tem_VS.jf > Tem_Vs.fa'
	if re.split('\.',read_1)[-1]=='gz':
		subprocess.check_output(cmd1,shell=True)
	else:
		os.system(cmd1)
	os.system(cmd2)
else:
	if re.split('\.',read_1)[-1]=='gz' or re.split('\.',read_2)[-1]=='gz':
		cmd1='zcat '+read_1+' '+read_2+' | '+file_dir+'/jellyfish-linux count /dev/fd/0 -m 25 -s 100M -t 8 --if '+snp_kmr_fa+' -o Tem_VS.jf '
	else:
		cmd1=file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if '+snp_kmr_fa+' -o Tem_VS.jf '+read_1+' '+read_2
	cmd2=file_dir+'/jellyfish-linux dump -c Tem_VS.jf > Tem_Vs.fa'
	if re.split('\.',read_1)[-1]=='gz' or re.split('\.',read_2)[-1]=='gz':
		subprocess.check_output(cmd1,shell=True)
	else:
		os.system(cmd1)
	os.system(cmd2)
	

freq_arr=[]
fnew=open('Tem_Vs.fa','r')
while True:
	line=fnew.readline().strip()
	if not line:break
	ele=line.split()
	dpsc[dkps[ele[0]]]+=int(ele[1])

carr=[]
for p in pos_snp:
	c=re.split('-',p)[0]
	if c not in carr:
		carr.append(c)
	if p not in dpsc:
		freq_arr.append(0)
	else:
		freq_arr.append(dpsc[p])

freq_arr=np.array(freq_arr)
#print(freq_arr)
#exit()

#freq_arr[freq_arr<=min_depth_absolute]=0
#keep=(freq_arr!=0)


### Check avg depth from the pos-snp frequency array
keep=(freq_arr!=0)
check_arr=freq_arr[keep]
if len(check_arr)==0:
	print('No kmers matched! No virus strain can be detected!')
	exit()
min_depth,max_depth=np.percentile(check_arr,[min_depth_percentile,max_depth_percentile])
#print(min_depth)
keep=np.logical_and.reduce((check_arr>=min_depth,check_arr<=max_depth))
check_arr2=check_arr[keep]
# Average depth of the frequency vector
min_depth_adf=min_depth_rate*np.mean(check_arr2)
##### !! freq arr filter using 2 or 0.1*avg-depth

if min_depth_adf<2:
	min_depth_adf=2 # Use 2 to test firstly
#min_depth_adf=min_depth_absolute
freq_arr[freq_arr<=min_depth_adf]=0
weighted_freq_arr=freq_arr/np.sum(freq_arr)

##### Filter done. ########


pos_freq_map=dict(zip(pos_snp,freq_arr))
#print(pos_freq_map['18722-C'])
t=np.array(pos_snp)
#sindex=np.argwhere(t=='18722-C')[0][0]
#print(freq_arr[sindex])
#exit()
#print(np.where(np.array(pos_snp)=='8946-T'))
	#ax=sns.distplot(a,norm_hist=False,kde=False)
	#ax=sns.distplot(a,norm_hist=False,kde=False)
#exit()

#print(pos_freq_map['8946-T'])
#exit()
#print(pos_freq_map['4248-T'])
#exit()
#des={}  # Strain -> the sum of frequency vector of this strain
ds_pos={}  # Strain -> 0 1 1 0 0 1 (pos-snp: yes or no vector)
ds_freq={} # Strain -> frequency vector of this strain
dmap_rate={} # Strain -> pos-snp map scpre of this strain
#ds_avgd={}  # Strain -> average depth of this strain
ds_num_m={} # Strain -> pos-snp map number
ds_num={} # Strain -> pos-snp map number, raw number
dmr={}  # Strain -> pos-snp map rate os this strain (map_number/raw_number)
#These 2 dict will be used to visualize the pos depth figure
dscf={} #Strain-> 110-A:110, 111-NA:0, ..
dscl={} # Strain-> 110-A:3000, 111-NA:0, ..
all_ps=[] # Record the matrix array
while True:
	line=f3.readline().strip()
	if not line:break
	ele=line.split('\t')
	#des[ele[0]]=0
	tem=[]
	#all_ps.append(ele[1:])
	for e in ele[1:]:
		tem.append(int(e))
	tem=np.array(tem)
	all_ps.append(tem)
	dscf[ele[0]]={}
	dscl[ele[0]]={}
	'''
	keep=(tem==1)
	ps_number=len(tem[keep])
	ds_num[ele[0]]=ps_number
	'''
	nt=freq_arr*tem # the frequency vector of this strain
	raw_c=len(tem[tem==1])
	map_c=len(nt[nt>0])
	'''
	if ele[0]=='>MT419847.1':
		print('>MT419847.1',map_c,raw_c)
	if ele[0]=='>MN938384.1':
		print('>MN938384.1',map_c,raw_c)
	'''
	map_rate=np.sum(tem*weighted_freq_arr)
	#map_rate=float(map_c)/float(raw_c)
	#dmap_rate[ele[0]]=map_rate
	dmap_rate[ele[0]]=map_rate
	ds_num[ele[0]]=str(map_c)+'/'+str(raw_c)
	ds_num_m[ele[0]]=int(map_c)
	if float(raw_c)==0:
		dmr[ele[0]]=0
	else:
		dmr[ele[0]]=float(map_c)/float(raw_c)
	#value=nt.sum()
	#des[ele[0]]=value
	ds_pos[ele[0]]=tem
	ds_freq[ele[0]]=nt
	## Get average depth of this strain ##
	'''
	keep=(nt!=0)
	nt=nt[keep]
	min_depth,max_depth=np.percentile(nt,[min_depth_percentile,max_depth_percentile])
	keep=np.logical_and.reduce((nt>=min_depth,nt<=max_depth))
	nt=nt[keep]
	ds_avgd[ele[0]]=np.mean(nt)
	'''
''' Here we initialize 2 dict for later visualization. '''
all_ps=np.array(all_ps)
all_sum=np.sum(all_ps,axis=0)
pos_label=dict(zip(pos_snp,list(all_sum)))
'''
for s in ds_freq:
	i=0
	for c in ds_freq[s]:
		if c==0:
			i+=1
			continue
		current_ps=pos_snp[i]	
		column=int(re.split('-',current_ps)[0])
		if column not in dscf[s]:
			dscf[s][column]=0
		if not c==0:
			dscf[s][column]+=c
		i+=1
	i2=0
	for c in ds_pos[s]:
		if c==0:
			i2+=1
			continue
		current_ps=pos_snp[i2]
		column=int(re.split('-',current_ps)[0])
		if column not in dscl[s]:
			dscl[s][column]=0
		if c==1:
			dscl[s][column]=pos_label[current_ps]
		i2+=1
print(len(dscl['>gb:J02176']),len(dscf['>gb:J02176']))
exit()
'''
#print(max(the_sum),len(np.argwhere(the_sum==len(all_ps_vec)-1)),len(np.argwhere(the_sum==1)))

#exit()
#index=np.where(np.array(pos_snp)=='8946-T')[0][0]
#print(ds_pos['>MT419847.1'][index])
#print(ds_pos['>MN938384.1'][index])
max_map=sorted(dmr.items(),key=lambda d:d[1],reverse=True)[0][1]
#exit()
if rks==1:
	res=sorted(ds_num_m.items(),key=lambda d:d[1],reverse=True)
else:
	res=sorted(dmap_rate.items(),key=lambda d:d[1],reverse=True)
top10_score_s=res[:10]
if rks==1:
	top10_score_s_old=top10_score_s
	top10_score_s=[]
	for r in top10_score_s_old:
		tem=(r[0],dmap_rate[r[0]])
		top10_score_s.append(tem)
#print(top10_score_s_old,top10_score_s)
#exit()

#print(top10_score_s)
#exit()
top_map_strain=[]
for r in res:
	if r[1]==res[0][1]:
		top_map_strain.append(r[0])
	else:break
#####  Unique nodes scan for all top strains
snp_arr=[]
### Pre-calculate the possible strain number, then decide whether should calculate weighted score.
pre_freq_arr=[]
#strain_num={}
for s in top_map_strain:
	snp_arr.append(ds_pos[s])
	pre_pa=ds_pos[s]*(-1)
	pre_pa=np.array(pre_pa)
	pre_pa[pre_pa==0]=1
	pre_freq_arr=freq_arr*pre_pa
	pre_freq_arr[pre_freq_arr<0]=0

keep=(pre_freq_arr!=0)
pre_freq_arr=pre_freq_arr[keep]

pre_pos_snp=np.array(pos_snp)[keep]
pre_ds_pos={}
for s in ds_pos:
	pre_ds_pos[s]=ds_pos[s][keep]
pre_wf_arr=pre_freq_arr/np.sum(pre_freq_arr)
strain_num={}
sn=0
#print(pre_freq_arr)
#exit()
while True:
	if len(pre_freq_arr)==0:break
	smr={}
	for r in ds_pos:
		if r in top_map_strain:continue
		tt=pre_ds_pos[r]
		if np.sum(tt)==0:continue
		#print(tt,pre_wf_arr)
		#exit()
		nt=tt*pre_wf_arr
		#mc=len(nt[nt>0])
		mr=np.sum(nt)
		smr[r]=mr
	res=sorted(smr.items(),key=lambda d:d[1],reverse=True)
	#print(res)
	#exit()
	ts=[]
	for r in res:
		if r[1]==res[0][1]:
			ts.append(r[0])
	#print(ts)
	if len(ts)==0:break
	if len(ts)>1:
		rmr={}
		for s in ts:
			rmr[s]=dmap_rate[s]
		res2=sorted(rmr.items(),key=lambda d:d[1],reverse=True)
		#for r in res2:
		strain_num[res2[0][0]]=''
	else:
		strain_num[ts[0]]=''
	#exit()
	#strain_num[res[0][0]]=''
	vm1=len(pre_freq_arr)
	pre_pa=pre_ds_pos[ts[0]]*(-1)
	pre_pa[pre_pa==0]=1
	pre_freq_arr=pre_freq_arr*pre_pa
	pre_freq_arr[pre_freq_arr<0]=0
	keep=(pre_freq_arr!=0)
	pre_freq_arr=pre_freq_arr[keep]
	if not np.sum(pre_freq_arr)==0:
		pre_wf_arr=pre_freq_arr/np.sum(pre_freq_arr)
	pre_pos_snp=pre_pos_snp[keep]
	for s in pre_ds_pos:
		pre_ds_pos[s]=pre_ds_pos[s][keep]
	vm=vm1-len(pre_freq_arr)
	if vm>1:
		#strain_num[res[0][0]]=''
		sn+=1
# Will recalculate the score and select top strain
if sn>1:
	#nw=len(strain_num)
	for s in top_map_strain:
		strain_num[s]=''
	sscore={}
	sna=[]
	for s in strain_num:
		sna.append(ds_pos[s])
	sna=np.array(sna)
	ssum=sna.sum(axis=0)
	ssum[ssum==0]=1
	for s in strain_num:
		snt=ds_pos[s]/ssum
		ns=dmap_rate[s]*snt
		ns=ns.sum(axis=0)
		sscore[s]=ns
	res=sorted(sscore.items(),key=lambda d:d[1],reverse=True)
	tem_map_strain=[]
	for r in res:
		if not dmr[r[0]]==max_map:continue
		tem_map_strain.append(r[0])
		break
	if len(tem_map_strain)>0:
		top_map_strain=tem_map_strain
	#print(sscore)

	
snp_arr=np.array(snp_arr)
pos_sum=snp_arr.sum(axis=0)
#pos_sum[pos_sum>1]=0
i=0

strain_unique={}
strain_unique_count={}
for p in pos_sum:
	column=pos_snp[i]
	#i+=1
	if p>=1:
		if pos_freq_map[column]<=min_depth_absolute:
			i+=1
			continue
		i2=0
		window=snp_arr[:,i]
		for w in window:
			if w>=1:
				if i2+1>len(top_map_strain):continue
				strain=top_map_strain[i2]
				if strain not in strain_unique:
					strain_unique[strain]={column:pos_freq_map[column]}
					strain_unique_count[strain]=1
				else:
					strain_unique[strain][column]=pos_freq_map[column]
					strain_unique_count[strain]+=1
			i2+=1
	i+=1
#print(freq_arr[1397])
#print(pos_freq_map['4248-T'])
#print(strain_unique)
#exit()
## Final output generate part ####
mp_strain=[]  # Most possible strain
op_strain=[]  # Other possible strain -> [S1,S2,S3,S5]
op_strain_batch=[] # Other possible strain ->[[S1,S2],[S3,S5],...]
#op_pos_snp=[] # Other possible pos-snp
#op_ps_strain=[] # The top map rate strain of other possible pos-snp
if not len(strain_unique)==0:
	#print(strain_unique)
	for s in strain_unique:
		mp_strain.append(s)
		#print(s,ds_avgd[s])
		#exit()
else:
	### Need to check whether these strains have close depth
	#d={}
	#strain_depth=[]
	for r in top_map_strain:
		mp_strain.append(r)
		#d[r]=ds_avgd[r]
		#raw_depth=ds_freq[r] # the freq arr of this strain -> [0, 22, 0 ...]
		#keep=(raw_depth>0)
		#nz_depth=raw_depth[keep]
		#print(np.mean(nz_depth))
		#exit()
	#exit()
	'''
	res2=sorted(d.items(),key=lambda d:d[1],reverse=True)
	for r in res2:
		mp_strain.append(r[0])
	'''
#print(strain_unique)
#exit()
#print(top_map_strain)
#for t in top_map_strain:
#print(t,dmap_rate[t])
## Check nodes of MT419847
## Iterative function to get other possible strains and other possible pos-snp
#print('Raw freq arr: ', freq_arr[sindex])
ds_avgd={}
for m in mp_strain:
	'''
	keep=(ds_freq[m]!=0)
	check_arr=ds_freq[m][keep]
	min_depth,max_depth=np.percentile(check_arr,[min_depth_percentile,max_depth_percentile])
	kp1=(check_arr<=min_depth)
	check_arr=check_arr[kp1]
	continue
	print(m,ds_freq[m][sindex])
	'''
	# Get the avg depth of the strain
	keep=(ds_freq[m]!=0)
	min_depth,max_depth=np.percentile(ds_freq[m][keep],[min_depth_percentile,max_depth_percentile])
	keep=np.logical_and.reduce((ds_freq[m]>=min_depth,ds_freq[m]<=max_depth))
	ds_avgd[m]=np.mean(ds_freq[m][keep])

	pos_arr=ds_pos[m]*(-1)
	pos_arr=np.array(pos_arr)
	pos_arr[pos_arr==0]=1
	freq_arr=freq_arr*pos_arr
	freq_arr[freq_arr<0]=0

keep=(freq_arr!=0)

left_freq_arr=freq_arr[keep]
#print('Filtered freq arr: ',left_freq_arr[sindex])
#exit()
pos_snp=np.array(pos_snp)
left_pos_snp=pos_snp[keep]
left_ds_pos={}
for s in ds_pos:
	left_ds_pos[s]=ds_pos[s][keep]
#print(left_freq_arr)
#print(left_pos_snp)
#exit()
#run=0
left_ps_freq_map=dict(zip(left_pos_snp,left_freq_arr))
#print(left_ps_freq_map)
#exit()
resl=sorted(left_ps_freq_map.items(),key=lambda d:d[1],reverse=True)
os_strain={} # '221-A-100-10000':['>MT312312.1',....]
os_arr=[] # ['221-A-100-10000','225-G-100-10000',....]
left_weighted_freq_arr=left_freq_arr/np.sum(left_freq_arr)
#print(left_ps_freq_map)
#print(left_freq_arr,left_pos_snp)
#exit()
#run=0
vmap={} # The dict used to record the valid map rate
#### Start Iterative process ######
if not len(left_freq_arr)==0:
	max_iter_times=len(left_freq_arr)
	for l in range(max_iter_times):
		#temd=dict(zip(left_pos_snp,left_freq_arr))
		#res=sorted(temd.items(),key=lambda d:d[1],reverse=True)
		#for r in res:
		if len(left_freq_arr)==0:break
		#column_snp=res[0][0]
		strain_map_rate={}
		for r in ds_pos:
			if r in mp_strain:continue
			#if r in op_strain:continue
			#if r in op_ps_strain:continue
			tem=left_ds_pos[r]
			nt=tem*left_weighted_freq_arr
			#raw_c=len(tem[tem==1])
			map_c=len(nt[nt>0])
			#map_rate=map_c
			map_rate=np.sum(nt)
			strain_map_rate[r]=map_rate
		selected_strain=[] # used to save the selected strains in this round
		res=sorted(strain_map_rate.items(),key=lambda d:d[1],reverse=True)
		#print(res)
		#exit()
		top_s=[]
		snp_arr=[]

		for r in res:
			if r[1]==res[0][1]:
				top_s.append(r[0])
				### Calculate avg depth
				nt=left_ds_pos[r[0]]*left_freq_arr
				keep=(nt!=0)
				nt=nt[keep]
				min_depth,max_depth=np.percentile(nt,[min_depth_percentile,max_depth_percentile])
				keep=np.logical_and.reduce((nt>=min_depth,nt<=max_depth))
				if not len(nt[keep])==0:
					nt=nt[keep]
				ds_avgd[r[0]]=np.mean(nt)
				## Done
				snp_arr.append(left_ds_pos[r[0]])
		#print(top_s)
		#exit()
		if len(top_s)>1:
			rank_map_rate={}
			for s in top_s:
				rank_map_rate[s]=dmap_rate[s]
			res=sorted(rank_map_rate.items(),key=lambda d:d[1],reverse=True)
			top_s=[]
			for r in res:
				top_s.append(r[0])
		top_s=np.array(top_s)
		pre=[]
		for r in resl:
			pre.append(r[0]+':'+str(r[1]))
		pre=','.join(pre)
		top_pos_snp=pre+'\t'+str(len(top_s))
		os_arr.append(top_pos_snp)
		os_strain[top_pos_snp]=[]
		for s in top_s:
			os_strain[top_pos_snp].append(s)
			pos_arr=left_ds_pos[s]*(-1)
			pos_arr[pos_arr==0]=1
			left_freq_arr=left_freq_arr*pos_arr
			left_freq_arr[left_freq_arr<0]=0
		keep=(left_freq_arr!=0)
		valid_map=len(left_freq_arr)
		left_freq_arr=left_freq_arr[keep]
		valid_map=valid_map-len(left_freq_arr)
		vmap[top_pos_snp]=valid_map
		if not np.sum(left_freq_arr)==0:
			left_weighted_freq_arr=left_freq_arr/np.sum(left_freq_arr)
		left_pos_snp=left_pos_snp[keep]
		left_ps_freq_map=dict(zip(left_pos_snp,left_freq_arr))
		for s in left_ds_pos:
			left_ds_pos[s]=left_ds_pos[s][keep]
		resl=sorted(left_ps_freq_map.items(),key=lambda d:d[1],reverse=True)

# All strain cross validation

### Output part
dc={}
if os.path.exists(db_dir+'/head_file.txt'):
	fp=open(db_dir+'/head_file.txt','r')
	dc={}
	while True:
		line=fp.readline().strip()
		if not line:break
		ele=line.split('\t')
		#pre=re.split('\|',line)[0].strip()
		#anno=re.split('\|',line)[1].strip()
		dc[ele[0]]=ele[1]

o=open(out_dir,'w+')
o.write('\t\tStrain_ID\tCls_info\tSubCls_info\tMap_Score\tValid_Map_Rate\tTotal_Map_Rate\tStrain_Depth\tStrain_info\tUnique_SNP\n')
o.write('>>Most possible strains:\n')
all_s=[]
vs_sd=[]
vs_so=[]
for s in mp_strain:
	all_s.append(s)
	vs_sd.append(s)
if len(os_arr)>0:
	for s in os_arr:
		all_s.append(os_strain[s][0])
		if vmap[s]>1:
			vs_so.append(os_strain[s][0])
# Strain-level identification 
s2cls={} # Used for output cluster info
s2sub={} # Used for output sub-cluster info
candidate_cls={}
fcls=open(cls_file)
while True:
	line=fcls.readline().strip()
	if not line:break
	ele=line.split('\t')
	if ele[0] not in all_s:
		s2cls[ele[0]]=ele[2]
		continue
	if ele[2]==ele[3]:
		s2cls[ele[0]]=ele[2]
		s2sub[ele[0]]='NA'
	else:
		s2cls[ele[0]]=ele[2]
		candidate_cls[ele[2]]=''
		#s2sub[ele[0]]=''
if not len(candidate_cls)==0:
	fsk=open(sub_kmr_file,'r')
	ksub={}
	cls_sub={}
	#subcount={}
	while True:
		line=fsk.readline().strip()
		if not line:break
		ele=line.split('\t')
		kmr=ele[0]
		if ele[1] in candidate_cls:
			if kmr not in ksub:
				ksub[kmr]={ele[1]:{}}
			if ele[1] not in ksub[kmr]:
				ksub[kmr][ele[1]]={}
			if ele[1] not in cls_sub:
				cls_sub[ele[1]]={}
			sub=re.split(',',ele[-1])
			for s in sub:
				ksub[kmr][ele[1]][s]=''
				cls_sub[ele[1]][s]=0
	ok=open('Tem_Vs2Sub.fa','w+')
	co=1
	for kmr in ksub:
		ok.write('>'+str(co)+'\n')
		ok.write(kmr+'\n')
		co+=1
	ok.close()
	if read_2=='':
		if re.split('\.',read_1)[-1]=='gz':
			cmd1='zcat '+read_1+' | '+file_dir+'/jellyfish-linux count /dev/fd/0 -m 25 -s 100M -t 8 --if Tem_Vs2Sub.fa -o Tem_VS2.jf'
		else:
			cmd1=file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if Tem_Vs2Sub.fa -o Tem_VS2.jf '+read_1
	else:
		if re.split('\.',read_1)[-1]=='gz' or re.split('\.',read_2)[-1]=='gz':
			cmd1='zcat '+read_1+' '+read_2+' | '+file_dir+'/jellyfish-linux count /dev/fd/0 -m 25 -s 100M -t 8 --if Tem_Vs2Sub.fa -o Tem_VS2.jf'
		else:
			cmd1=file_dir+'/jellyfish-linux count -m 25 -s 100M -t 8 --if Tem_Vs2Sub.fa -o Tem_VS2.jf '+read_1+' '+read_2
	cmd2=file_dir+'/jellyfish-linux dump -c Tem_VS2.jf > Tem_Vs2.fa'
	if re.split('\.',read_1)[-1]=='gz':
		subprocess.check_output(cmd1,shell=True)
	else:
		os.system(cmd1)
	os.system(cmd2)
	ft2=open('Tem_Vs2.fa','r')
	while True:
		line=ft2.readline().strip()
		if not line:break
		ele=line.split()
		if int(ele[1])>=min_depth_adf:
			for c in ksub[ele[0]]:
				for c2 in ksub[ele[0]][c]:
					cls_sub[c][c2]+=int(ele[1])
	for s in all_s:
		if s in s2sub:continue
		if s2cls[s] not in cls_sub:
			s2sub[s]='NA'
		else:
			res=sorted(cls_sub[s2cls[s]].items(),key=lambda d:d[1],reverse=True)
			if len(res)>1:
				if res[0][1]==res[1][1]:
					s2sub[s]='NA'
				else:
					s2sub[s]=res[0][0]
			else:
				s2sub[s]=res[0][0]


for s in mp_strain:
	#keep=(ds_pos[s]==1)
	#all_s.append(s)
	if s in strain_unique:
		if s not in dc:
			o.write('\t\t'+s+'\t'+s2cls[s]+'\t'+s2sub[s]+'\t'+str(dmap_rate[s])+'\t'+str(ds_num[s])+'\t'+str(ds_num[s])+'\t'+str(ds_avgd[s])+'\t\t\t'+str(strain_unique[s])+'\n')
		else:
			o.write('\t\t'+s+'\t'+s2cls[s]+'\t'+s2sub[s]+'\t'+str(dmap_rate[s])+'\t'+str(ds_num[s])+'\t'+str(ds_num[s])+'\t'+str(ds_avgd[s])+'\t'+dc[s]+'\t'+str(strain_unique[s])+'\n')
	else:
		if s not in dc:
			o.write('\t\t'+s+'\t'+s2cls[s]+'\t'+s2sub[s]+'\t'+str(dmap_rate[s])+'\t'+str(ds_num[s])+'\t'+str(ds_num[s])+'\t'+str(ds_avgd[s])+'\t\t\tNA\n')
		else:
			o.write('\t\t'+s+'\t'+s2cls[s]+'\t'+s2sub[s]+'\t'+str(dmap_rate[s])+'\t'+str(ds_num[s])+'\t'+str(ds_num[s])+'\t'+str(ds_avgd[s])+'\t'+dc[s]+'\tNA\n')
o.write('>>Other possible strains:\n')
if len(os_arr)>0:
	#i=1
	for s in os_arr:
		ele=re.split('\t',s)
		if vmap[s]==1:
			o.write('\t>>(Could be FP)'+ele[0]+',Genome_num: '+ele[1]+'\n')
		else:
			o.write('\t>>'+ele[0]+',Genome_num: '+ele[1]+'\n')
		all_s.append(os_strain[s][0])
		for n in os_strain[s]:
			#all_s.append(n)
			a=re.split('/',ds_num[n])[-1]
			vm=str(vmap[s])+'/'+a
			if n in s2sub:
				if n not in dc:
					o.write('\t\t'+n+'\t'+s2cls[n]+'\t'+s2sub[n]+'\t'+str(dmap_rate[n])+'\t'+vm+'\t'+ds_num[n]+'\t'+str(ds_avgd[n])+'\t\t\tNot_record\n')
				else:
					o.write('\t\t'+n+'\t'+s2cls[n]+'\t'+s2sub[n]+'\t'+str(dmap_rate[n])+'\t'+vm+'\t'+ds_num[n]+'\t'+str(ds_avgd[n])+'\t'+dc[n]+'\tNot_record\n')
			else:
				if n not in dc:
					o.write('\t\t'+n+'\t'+s2cls[n]+'\tNot_record\t'+str(dmap_rate[n])+'\t'+vm+'\t'+ds_num[n]+'\t'+str(ds_avgd[n])+'\t\t\tNot_record\n')
				else:
					o.write('\t\t'+n+'\t'+s2cls[n]+'\tNot_record\t'+str(dmap_rate[n])+'\t'+vm+'\t'+ds_num[n]+'\t'+str(ds_avgd[n])+'\t'+dc[n]+'\tNot_record\n')
				
		#outs='\t'.join(os_strain[s])
		#o.write('\t\t\t'+outs+'\n')
		#i+=1
else:
	o.write('\tCan not detect other strains.\n')
res=sorted(dmr.items(),key=lambda d:d[1],reverse=True)
o.write('\n>>Highest_Map_Strains (Could be FP):\n')
final={}
for r in res:
	if r[1]==res[0][1]:
		if r[0] not in all_s:
			final[r[0]]=dmap_rate[r[0]]
if not len(final)==0:
	res=sorted(final.items(),key=lambda d:d[1],reverse=True)
	for s in res:
		if s[0] in s2sub:
			if s[0] not in dc:
				o.write('\t\t'+s[0]+'\t'+s2cls[s[0]]+'\t'+s2sub[s[0]]+'\t'+str(dmap_rate[s[0]])+'\t'+ds_num[s[0]]+'\t'+ds_num[s[0]]+'\tNA\tNA\n')
			else:
				o.write('\t\t'+s[0]+'\t'+s2cls[s[0]]+'\t'+s2sub[s[0]]+'\t'+str(dmap_rate[s[0]])+'\t'+ds_num[s[0]]+'\t'+ds_num[s[0]]+'\t'+dc[s[0]]+'\tNA\n')
		else:
			if s[0] not in dc:
				o.write('\t\t'+s[0]+'\t'+s2cls[s[0]]+'\tNot_record\t'+str(dmap_rate[s[0]])+'\t'+ds_num[s[0]]+'\t'+ds_num[s[0]]+'\tNA\tNA\n')
			else:
				o.write('\t\t'+s[0]+'\t'+s2cls[s[0]]+'\tNot_record\t'+str(dmap_rate[s[0]])+'\t'+ds_num[s[0]]+'\t'+ds_num[s[0]]+'\t'+dc[s[0]]+'\tNA\n')

o.write('>>Top10_Score_Strains:\n')
for t in top10_score_s:
	if t[0] in s2sub:
		if t[0] not in dc:
			o.write('\t\t'+t[0]+'\t'+s2cls[t[0]]+'\t'+s2sub[t[0]]+'\t'+str(t[1])+'\t'+ds_num[t[0]]+'\t'+ds_num[t[0]]+'\tNA\tNA\n')
		else:
			o.write('\t\t'+t[0]+'\t'+s2cls[t[0]]+'\t'+s2sub[t[0]]+'\t'+str(t[1])+'\t'+ds_num[t[0]]+'\t'+ds_num[t[0]]+'\t'+dc[t[0]]+'\tNA\n')
	else:
		if t[0] not in dc:
			o.write('\t\t'+t[0]+'\t'+s2cls[t[0]]+'\tNot_record\t'+str(t[1])+'\t'+ds_num[t[0]]+'\t'+ds_num[t[0]]+'\tNA\tNA\n')
		else:
			o.write('\t\t'+t[0]+'\t'+s2cls[t[0]]+'\tNot_record\t'+str(t[1])+'\t'+ds_num[t[0]]+'\t'+ds_num[t[0]]+'\t'+dc[t[0]]+'\tNA\n')
## Remove tem file 
os.system('rm Tem_Vs* Tem_VS*')
## From this line, we will generate strain-level analysis report
print('Txt report is done. Now will generate pdf report!')

vs_so=vs_so[:5]
for s in ds_freq:
	check=0
	if s  in vs_sd:
		check=1
	if s in vs_so:
		check=2
	if check==0:continue
	i=0
	for c in ds_freq[s]:
		if c==0:
			i+=1
			continue
		current_ps=pos_snp[i]
		column=int(re.split('-',current_ps)[0])
		if column not in dscf[s]:
			dscf[s][column]=c
		
		i+=1
	i2=0
	for c in ds_pos[s]:
		if c==0:
			i2+=1
			continue
		current_ps=pos_snp[i2]
		column=int(re.split('-',current_ps)[0])
		if column not in dscl[s]:
			dscl[s][column]=0
		if c==1:
			dscl[s][column]=pos_label[current_ps]
		i2+=1

ov1=open('Mps_ps_depth.csv','w+')
ov2=open('Ops_ps_depth.csv','w+')
ov1.write('ID,Column_ID')
#carr=[]
for s in vs_sd:
	ov1.write(','+s+'_Freq')
	ov1.write(','+s+'_LNum')
ov1.write('\n')
#carr=sorted(list(dscf[vs_sd[0]].keys()))
carr=[]
for c in pos_snp:
	c=re.split('-',c)[0]
	if int(c) not in carr:
		carr.append(int(c))

i=1
for c in carr:
	ov1.write(str(i)+','+str(c))
	check=0
	check1=0
	check2=0
	for s in vs_sd:
		if c not in dscf[s]:
			check+=1
			check1=1
		if c not in dscl[s]:
			check+=1
			check2=1
		if check==0:
			ov1.write(','+str(dscf[s][c])+','+str(dscl[s][c]))
		elif check==1:
			if check1==1:
				ov1.write(',0,'+str(dscl[s][c]))	
	ov1.write('\n')
	i+=1
	
	

ov2.write('ID,Column_ID')
if len(vs_so)==0:
	ov2.write(',None,None\n')
else:
	for s in vs_so:
		ov2.write(','+s+'_Freq')
		ov2.write(','+s+'_LNum')
	ov2.write('\n')
	i=1
	for c in carr:
		ov2.write(str(i)+','+str(c))
		check=0
		check1=0
		check2=0
		for s in vs_so:
			if s not in dscf or s not in dscl:
				print('Warning: ',s,' not in final dict!')
				continue
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
