import re
import os
import getopt
import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

input_msa=''
k=25  # default kmer size
cfile='' # Column file generated by perl file from Prof Sun and transformed by my code
out='Pos-snp-kmer-all.txt' # default
out2='Pos-snp-kmer-all.fa' # default

opts,args=getopt.getopt(sys.argv[1:],"i:k:c:o:f:")

for opt,arg in opts:
	if opt=='-i':
		input_msa=arg
	elif opt=='-k':
		k=int(arg)
	elif opt=='-c':
		cfile=arg
	elif opt=='-o':
		out=arg
	elif opt=='-f':
		out2=arg


dmap={'a':'A','t':'T','g':'G','c':'C','A':'A','T':'T','G':'G','C':'C'}
check_dict={'a':'A','t':'T','g':'G','c':'C','A':'A','T':'T','G':'G','C':'C','-':''}

f1=open(cfile,'r')
dseq={}
column_arr=[]
#dpbk={} #dict -> pos_base_kmr
while True:
	line=f1.readline().strip()
	if not line:break
	if not re.search('column',line):continue
	ele=line.split()
	'''
	if len(ele)<3:
		print(ele)
		exit()
	'''
	if int(ele[1]) not in column_arr:
		column_arr.append(int(ele[1]))
	'''
	#if ele[-2] not in dpbk:
		#dpbk[ele[-2]]={dmap[ele[-1]]:{}}
	#else:
		#if dmap[ele[-1]] in dpbk[ele[-2]]:
			#print('Double Column SNP - Check plz.')
			#exit()
		#dpbk[ele[-2]][dmap[ele[-1]]]={}
	for e in ele[1:-2]:
		if re.search(':',e):
			name=re.split(':',e)[0]
		else:
			name=e
		dseq[name]=''
	'''
#print(len(column_arr))
column_arr=sorted(column_arr)

f2=open(input_msa,'r')
dmsa={}
ds2id={}   # Strain_Name -> Strain_ID
did2s={}
oi=open('ID2Name.txt','w+')
i=0
while True:
	line=f2.readline().strip()
	if not line:break
	if re.search('>',line):
		name=line.split()[0].strip()
		#if name not in dseq:continue
		dmsa[name]=''
		oi.write(str(i)+'\t'+name+'\n')
		ds2id[name]=str(i)
		did2s[i]=name
		i+=1
	else:
		#if name not in dseq:continue
		dmsa[name]+=line
'''
ot=open('Target_msa.fasta','w+')
target=[3344,0,9547,9548,9549,9550]
for t in target:
	ot.write('>'+did2s[t]+'\n')
	s=''
	for c in column_arr:
		s+=dmsa[did2s[t]][c]
	ot.write(s+'\n')
exit()
'''
### Most important part ###
dpbk={} # Dict -> Pos:Base:Kmr
dkpb={}	# Dict -> Kmr: POS-Base
dkl={}  # Dict -> Kmr : Label
dks={}  # Dict -> kmr : Strain
total_column=len(column_arr)
i=1
for c in column_arr:
	print('Progress: ',i,'/',total_column)
	i+=1
	sid=0
	for s in dmsa:
		if dmsa[s][c] not in dmap:continue
		base=dmap[dmsa[s][c]]
		left=dmsa[s][:c][::-1]  # Left Sequence of column base
		right=dmsa[s][c+1:]	    # Right Sequence of column base
		pos_base=str(c)+'-'+base
		block_seq=dmsa[s][c-12:c]+base+dmsa[s][c+1:c+13]
		check_num=0
		for b in block_seq:
			if b not in check_dict:
				check_num=1
				break
		if check_num==1:
			continue

		#print(len(dmsa[s][c-12:c]),len(dmsa[s][c+1:c+13]))
		#exit()

		lseq=''
		rseq=''
		for l in left:
			if len(lseq)==(k-1)/2:
				break
			if l in dmap:
				lseq+=dmap[l]
			elif not l=='-':
				lseq+=l.upper()
		for r in right:
			if len(rseq)==(k-1)/2:
				break
			if r in dmap:
				rseq+=dmap[r]
			elif not r=='-':
				rseq+=r.upper()
		lseq=lseq[::-1]
		kmr=lseq+base+rseq
		if len(kmr) <k:
			if len(lseq)<(k-1)/2:
				rseq=''
				rl=k-1-len(lseq)
				for r in right:
					if len(rseq)==rl:
						break
					if r in dmap:
						rseq+=dmap[r]
					elif not r=='-':
						rseq+=r.upper()
			if len(rseq)<(k-1)/2:
				lseq=''
				ll=k-1-len(rseq)
				for l in left:
					if len(lseq)==ll:
						break
					if l in dmap:
						lseq+=dmap[l]
					elif not l=='-':
						lseq+=l.upper()
			kmr=lseq+base+rseq
		if c not in dpbk:
			dpbk[c]={base:{kmr:1}}
		else:
			if base not in dpbk[c]:
				dpbk[c][base]={kmr:1}
			else:
				dpbk[c][base][kmr]=1
		if kmr not in dkl:
			dkl[kmr]={ds2id[s]:''}
		else:
			dkl[kmr][ds2id[s]]=''
		if kmr not in dkpb:
			dkpb[kmr]={pos_base:''}
		else:
			dkpb[kmr][pos_base]=''

o=open(out,'w+')
o3=open(out2,'w+')
i=1
for k in dkpb:
	info=','.join(list(dkpb[k].keys()))
	info2=','.join(list(dkl[k].keys()))
	o.write(k+'\t'+info+'\t'+info2+'\n')
	o3.write('>'+str(i)+'\n'+k+'\n')
	seq=Seq(k,IUPAC.ambiguous_dna)
	i+=1
	k2=seq.reverse_complement()._data
	o.write(k2+'\t'+info+'\t'+info2+'\n')
	o3.write('>'+str(i)+'\n'+k2+'\n')
	i+=1
print('Kmer POS-Snp file done! Start generating matrix.')
#### Generate Strain pos-snp matrix ####

o2=open('Strain_pos_snp_matrix_consider_all.txt','w+')
base_arr=['A','T','G','C']
pos_snp=[]
for c in column_arr:
	for b in base_arr:
		pb=str(c)+'-'+b
		pos_snp.append(pb)
fl='\t'.join(pos_snp)
o2.write('\t\t'+fl+'\n')
for s in dmsa:
	temd={}
	for c in column_arr:
		if dmsa[s][c] not in dmap:
			temd[str(c)+'-N']=''
		else:
			temd[str(c)+'-'+dmap[dmsa[s][c]]]=''
	o2.write(s)
	for p in pos_snp:
		if p not in temd:
			o2.write('\t0')	
		else:
			o2.write('\t1')
	o2.write('\n')

