import re
import os
import getopt
import sys
from Bio import AlignIO
import collections
import math


input_msa=''
dash_cutoff=0.01
opts,args=getopt.getopt(sys.argv[1:],"i:")

for opt,arg in opts:
	if opt=='-i':
		input_msa=arg

def estimate_shannon_entropy(dna_seq):
	#m=len(dna_seq)
	m=0
	base_dict={'a':'','t':'','g':'','c':'','A':'','T':'','G':'','C':''}
	bases_raw=dict(collections.Counter([tmp_base for tmp_base in dna_seq]))
	''' we will not consider dash or special characters when we calculate entropy '''
	bases={}
	for b in bases_raw:
		if b not in base_dict:continue
		bases[b]=bases_raw[b]
		m+=bases_raw[b]
	shannon_entropy_value=0
	for base in bases:
		n_i=bases[base]
		p_i=n_i / float(m)
		entropy_i=p_i * (math.log(p_i))
		shannon_entropy_value+=entropy_i
	return shannon_entropy_value * (-1), bases_raw
def estimate_shannon_entropy_iterate(alignment,used_strains,used_columns,centropy):
	new_centropy={}
	new_cbase_freq={}
	base_dict={'a':'','t':'','g':'','c':'','A':'','T':'','G':'','C':''}
	for c in centropy:
		tem_bf={}
		m=0
		if c in used_columns:
			new_centropy[c]=1000
			continue
		else:
			i=0
			for s in alignment[:,c]:
				if i in used_strains:
					i+=1
					continue
				if s not in base_dict:
					i+=1
					continue
				if s not in tem_bf:
					tem_bf[s]=1
					m+=1
					i+=1
				else:
					tem_bf[s]+=1
					m+=1
					i+=1
		shannon_entropy_value=0
		for base in tem_bf:
			n_i=tem_bf[base]
			p_i=n_i / float (m)
			entropy_i=p_i * (math.log(p_i))
			shannon_entropy_value+=entropy_i
		new_centropy[c]=shannon_entropy_value * (-1)
		new_cbase_freq[c]=tem_bf
	return new_centropy, new_cbase_freq
 
		

def cls_with_entropy(alignment, centropy ,cbase_freq, sid, sname2seq, sid2name, sname2id):
	used_strains={}
	used_columns={}
	#copyentropy=centropy # copyentropy will be updated during the iterative process
	#copyalignment=alignment # copyalignment will not be updated
	base_dict={'a':'','t':'','g':'','c':'','A':'','T':'','G':'','C':''}
	o1=open('output_test.txt','w+')
	o1.write('finish inputting '+str(len(sname2id))+' sequences from the alignment\n')
	cqnum=len(centropy)
	snum=len(sname2id)
	while True:
		if len(used_strains)==len(sname2id) or len(used_columns)==cqnum:break
		print('Progress: C: ',len(used_columns),'/',cqnum,' S:',len(used_strains),'/',snum)
		res=sorted(centropy.items(),key=lambda d:d[1])
		## Select columns with minimum entropy
		for r in res:	
			if r[0] not in used_columns: # If the column is already used, we will not consider it again
				current_c=r[0]	# This is the column id, the column with minimum entropy
				break
		rawseq=alignment[:,current_c]
		#base_freq_raw=dict(collections.Counter([tmp_base for tmp_base in seq]))
		#base_freq={}
		'''
		for b in base_freq_raw:
			if b in base_dict:
				base_freq[b]=base_freq_raw[b]
		'''
		res2=sorted(cbase_freq[current_c].items(),key=lambda d:d[1])
		round_strain=[]	
		for r in res2:
			if r[0] not in base_dict:continue
			if r[0]==res2[-1][0]:continue
			strain_id=[i for i, letter in enumerate(rawseq) if letter == r[0]]
			round_strain=round_strain+strain_id
			# Output cls info
			o1.write('column '+str(current_c))
			for b in ['a','t','g','c']:
				if b not in cbase_freq[current_c]:
					o1.write(' '+b+'(0)')
				else:
					o1.write(' '+b+'('+str(cbase_freq[current_c][b])+')')
			o1.write('\n')
			for s in strain_id:
				o1.write('>'+sid2name[s]+' ')
			o1.write(str(current_c)+' '+r[0]+'\n')
			cbase_freq[current_c][r[0]]=0
		for s in round_strain:
			used_strains[s]=''
		# record used columns 
		used_columns[current_c]=''
		# update entropy and cbase_freq
		centropy, cbase_freq=estimate_shannon_entropy_iterate(alignment,used_strains,used_columns,centropy)
		
		
# Load msa
alignment=AlignIO.read(input_msa,'fasta')
#entropy,base_count=estimate_shannon_entropy(alignment[:,7])
sid=range(len(alignment[:,0]))
sname=[]
sname2seq={}  # Strain_Name -> Strain_msa_sequence
cnum=0
for record in alignment:
	name=record.id.strip()
	sname.append(name)	
	cnum=len(record.seq)
	sname2seq[name]=record.seq
cid=range(cnum)
centropy={} # Column_Entropy
cid_seq={} # Column_ID -> Column_Sequence
cbase_freq={}
# Here we can filter some columns
for column in cid:
	seq=alignment[:,column]
	entropy,base_count=estimate_shannon_entropy(seq)
	if '-' not in base_count:
		if not entropy==0:
			centropy[column]=entropy 
			#cid_seq[column]=seq
			cbase_freq[column]=base_count
	#res=sorted(raw_base.items(),key=lambda d:d[1],reverse=True)
	else:
		if base_count['-']<= dash_cutoff * len(sid):
			#centropy.append(entropy)
			if not entropy==0:
				#cid_seq[column]=seq
				cbase_freq[column]=base_count
				centropy[column]=entropy
sid2name=dict(zip(sid,sname))  # Strain_ID -> Strain_Name
sname2id=dict(zip(sname,sid))  # Strain_Name -> Strain_ID
### Now we need to iterate all columns
cls_with_entropy(alignment, centropy, cbase_freq, sid, sname2seq, sid2name, sname2id)
