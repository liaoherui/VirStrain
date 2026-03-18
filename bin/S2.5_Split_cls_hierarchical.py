import re
import os
import getopt
import sys
from Bio import AlignIO
import collections
import math
import matplotlib.pyplot as plt
import networkx
from networkx.algorithms.components.connected import connected_components
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


input_msa=''
column_file=''
cls_file=''
dash_cutoff=0.01
k=25
opts,args=getopt.getopt(sys.argv[1:],"i:c:r:")

for opt,arg in opts:
	if opt=='-i':
		input_msa=arg
	if opt=='-c':
		column_file=arg
	if opt=='-r':
		cls_file=arg	

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

def estimate_shannon_entropy_spe(alignment,strain_set,left_c):
	m=0
	base_dict={'a':'','t':'','g':'','c':'','A':'','T':'','G':'','C':''}
	centropy={}
	cbase_freq={}
	for c in left_c:
		seq=''
		i=0
		for n in alignment[:,c]:
			if i not in strain_set:
				i+=1
				continue
			seq+=n
			i+=1
		entropy,base_count=estimate_shannon_entropy(seq)
		cbase_freq[c]=base_count
		if entropy==0:
			centropy[c]=1000
		else:
			centropy[c]=entropy
	#print(centropy,cbase_freq)
	return centropy,cbase_freq
	
def estimate_shannon_entropy_iterate(alignment,used_strains,used_columns,centropy,strain_set):
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
				if i not in strain_set:
					i+=1
					continue
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
 
		

#def cls_with_entropy_hier(alignment, centropy ,cbase_freq, sid, sname2seq, sid2name, sname2id):
def cls_with_entropy_hier(alignment,left_c,strain_set):
	used_strains={}
	used_columns={}
	selected_columns={}
	#copyentropy=centropy # copyentropy will be updated during the iterative process
	#copyalignment=alignment # copyalignment will not be updated
	base_dict={'a':'','t':'','g':'','c':'','A':'','T':'','G':'','C':''}
	#o1=open('output_test.txt','w+')
	#o1.write('finish inputting '+str(len(sname2id))+' sequences from the alignment\n')
	cqnum=len(left_c)
	snum=len(strain_set)
	sarr=sorted(list(strain_set.keys()))
	centropy, cbase_freq=estimate_shannon_entropy_spe(alignment,strain_set,left_c)
	'''
	if snum==66:
		print(centropy)
		print(cbase_freq)
		exit()
	'''
	for c in centropy:
		if centropy[c]==1000:
			used_columns[c]=''
	#exit()
	ibreak=0
	while True:
		if len(used_strains)==snum or len(used_columns)==cqnum:break
		print('Progress: C: ',len(used_columns),'/',cqnum,' S:',len(used_strains),'/',snum)
		res=sorted(centropy.items(),key=lambda d:d[1])
		'''
		if snum==66:
			print(res)
		'''
		if res[0][1]==1000 and ibreak==0:
			return {},{},{},[]		
		ibreak+=1
		if res[0][1]==1000:break
		## Select columns with minimum entropy
		for r in res:	
			if r[0] not in used_columns: # If the column is already used, we will not consider it again
				current_c=r[0]	# This is the column id, the column with minimum entropy
				break
		rawseq=''
		sid2ni={}
		i=0
		for s in sarr:		
			rawseq+=alignment[s,current_c]
			sid2ni[i]=s
			i+=1
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
			strain_id=[sid2ni[i] for i, letter in enumerate(rawseq) if letter == r[0]]
			round_strain=round_strain+strain_id
			# Output cls info
			if not centropy[current_c]==1000:
				selected_columns[current_c]=''
			'''
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
			'''
			cbase_freq[current_c][r[0]]=0
		for s in round_strain:
			used_strains[s]=''
		# record used columns 
		used_columns[current_c]=''
		# update entropy and cbase_freq
		centropy, cbase_freq=estimate_shannon_entropy_iterate(alignment,used_strains,used_columns,centropy,strain_set)
	## Cluster based on selected column
	sub_info={} # Sub-cluster-ID -> Strain_Name
	kmr_sub={} # kmer -> sub-cluster
	kmr_pos={} # kmer -> pos-snp
	strain_sub={}
	sort_sub=[]
	if len(selected_columns)==0:
		return 	sub_info,kmr_sub,kmr_pos,sort_sub
	else:
		# For cluster
		def to_graph(l):
			G=networkx.Graph()
			for part in l:
				G.add_nodes_from(part)
				G.add_edges_from(to_edges(part))
			return G
		def to_edges(l):
			it=iter(l)
			last=next(it)
			for current in it:
				yield last,current
				last=current
	
		carr=sorted(list(selected_columns.keys()))
		scseq=[]
		for s in strain_set:
			tem=[]
			cseq=''
			for c in carr:
				cseq+=alignment[s,c]
			tem.append(cseq)
			tem_s=str(s)+'_cls'
			tem.append(tem_s)
			scseq.append(tem)
		#print(scseq)
		#exit()
		G=to_graph(scseq)
		count=1
		for ele in list(connected_components(G)):
			ele=list(ele)
			pre='SubCls'+str(count)+'_'+str(len(ele)-1)
			sort_sub.append(pre)
			sub_info[pre]=[]
			for e in ele:
				if not re.search('_',e):continue
				e=re.sub('_cls','',e)
				e=int(e)
				sub_info[pre].append(e)
				strain_sub[e]=pre
			count+=1
		# Ready to extract kmer for each sub-cluster
		check_dict={'a':'A','t':'T','g':'G','c':'C','A':'A','T':'T','G':'G','C':'C','-':''}
		dmap={'a':'A','t':'T','g':'G','c':'C','A':'A','T':'T','G':'G','C':'C'}
		for c in carr:
			for s in strain_set:
				base=alignment[s,c].upper()
				left=alignment[s].seq[:c][::-1]
				right=alignment[s].seq[c+1:]
				pos_base=str(c)+'-'+base
				if c-12<0:
					block_seq=alignment[s].seq[0:c]+base+alignment[s].seq[c+1:c+13]
				else:
					'''
					print(c+13,c-12,base)
					print(alignment[s].seq[c-12:c])
					print(alignment[s].seq[c+1:c+13])
					print(len(alignment[s,:]))
					exit()
					'''
					block_seq=alignment[s].seq[c-12:c]+base+alignment[s].seq[c+1:c+13]
				check_num=0
				for b in block_seq:
					if b not in check_dict:
						check_num=1
						break
				if check_num==1:continue
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
				if len(kmr) <25:
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
				if kmr not in kmr_sub:
					kmr_sub[kmr]={strain_sub[s]:''}
				else:
					kmr_sub[kmr][strain_sub[s]]=''
				if kmr not in kmr_pos:
					kmr_pos[kmr]={pos_base:''}
				else:
					kmr_pos[kmr][pos_base]=''
		# recluster according to pos-snp kmer
		#for k in kmr_sub:
		
		return  sub_info,kmr_sub,kmr_pos,sort_sub				
# Load used column info
cfile=open(column_file,'r')
bused_c={}
while True:
	line=cfile.readline().strip()
	if not line:break
	if not re.search('column',line):continue
	ele=line.split()
	bused_c[int(ele[1])]=''
# Load msa
alignment=AlignIO.read(input_msa,'fasta')
#entropy,base_count=estimate_shannon_entropy(alignment[:,7])
sid=range(len(alignment[:,0]))
sname=[]
sname2seq={}  # Strain_Name -> Strain_msa_sequence
cnum=0
for record in alignment:
	name=record.id.strip()
	sname.append('>'+name)	
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
## Check how many columns left
stat=0
left_c={}
for c in centropy:
	if c not in bused_c:
		left_c[c]=''
		stat+=1
print('Total columns: ',cnum)
print('Used column in first iterate: ',len(bused_c))
print('Left Column: ',stat,'/',len(centropy))
## We need to parse the cluster, and exclude used strains
#finished_strains={} # ID -> Strain_Name
strain_info={} # Strain_Name -> Cluster_info
cls_strain={} # Cluster -> Strain_Name -> Strain_ID
cls_rep={}
strain_rep={} # Strain_Name -> Rep_Strain
fc=open(cls_file,'r')
while True:
	line=fc.readline().strip()
	if not line:break
	if re.search('Cluster',line):
		cls=re.sub('>','',line)
		cls_strain[cls]={}
	else:
		ele=line.split('\t')
		if re.search('\*',line):
			rep=ele[0]
			#cls_rep[cls]=rep
			strain_rep[rep]=rep
			strain_info[rep]=cls
			cls_strain[cls][rep]=int(ele[-1])
			cls_rep[cls]=rep
		else:
			strain_rep[ele[0]]=rep
			strain_info[ele[0]]=cls
			cls_strain[cls][ele[0]]=int(ele[-1])

#print(cls_strain['Cluster14_2'],strain_rep['>gb:CY067031'],strain_info['>gb:CY067031'])

#exit()
### Most important part ----> Split cluster
o1=open('Strain_cls_info.txt','w+')
o2=open('Rebuild_cls.clstr','w+')
o3=open('SubCls_kmer.txt','w+')
raw_cls_size=[]
new_cls_size=[]
for c in cls_strain:
	raw_cls_size.append(len(cls_strain[c]))
	if len(cls_strain[c])==1:
		new_cls_size.append(len(cls_strain[c]))
		for s in cls_strain[c]:
			o1.write(s+'\t'+strain_rep[s]+'\t'+strain_info[s]+'\t'+strain_info[s]+'\t'+str(sname2id[s])+'\t1\n')
			o2.write('>'+c+'\n')
			o2.write('\t'+s+'\t*\t'+str(sname2id[s])+'\n')
	else:
		strain_set={}
		for s in  cls_strain[c]:
			strain_set[cls_strain[c][s]]=''
		sub_info,kmr_sub,kmr_pos,sort_sub=cls_with_entropy_hier(alignment,left_c,strain_set)
		print(sub_info, sort_sub)
		#exit()
		if len(sub_info)==0 or len(sort_sub)==1:
			new_cls_size.append(len(cls_strain[c]))
			o2.write('>'+c+'\n')
			for s in cls_strain[c]:
				o1.write(s+'\t'+strain_rep[s]+'\t'+strain_info[s]+'\t'+strain_info[s]+'\t'+str(sname2id[s])+'\t'+str(len(cls_strain[c]))+'\n')
				if s==cls_rep[c]:
					o2.write('\t'+s+'\t*\t'+str(sname2id[s])+'\n')
				else:
					o2.write('\t'+s+'\t'+str(sname2id[s])+'\n')
		else:
			o2.write('>'+c+'\n')
			for sub in sort_sub:
				o2.write('\t>'+sub+'\n')
				new_cls_size.append(len(sub_info[sub]))
				for ss in sub_info[sub]:
					ss=sid2name[ss]
					o1.write(ss+'\t'+strain_rep[ss]+'\t'+strain_info[ss]+'\t'+sub+'\t'+str(sname2id[ss])+'\t'+str(len(sub_info[sub]))+'\n')
					if ss==cls_rep[c]:
						o2.write('\t\t'+ss+'\t*\t'+str(sname2id[ss])+'\n')
					else:
						o2.write('\t\t'+ss+'\t'+str(sname2id[ss])+'\n')
			for kl in kmr_sub:
				info1=','.join(list(kmr_pos[kl].keys()))
				info2=','.join(list(kmr_sub[kl].keys()))
				o3.write(kl+'\t'+c+'\t'+info1+'\t'+info2+'\n')
				seq=Seq(kl,IUPAC.ambiguous_dna)
				k2=seq.reverse_complement()._data
				o3.write(k2+'\t'+c+'\t'+info1+'\t'+info2+'\n')
				

plt.hist(raw_cls_size,bins=100)
plt.xlabel('Cluster size')
plt.ylabel('Cluster Number')
plt.savefig('Before_split.png')	
plt.figure()
plt.hist(new_cls_size,bins=100)
plt.xlabel('Cluster size')
plt.ylabel('Cluster Number')
plt.savefig('After_split.png')
### Now we need to iterate all columns
#cls_with_entropy(alignment, centropy, cbase_freq, sid, sname2seq, sid2name, sname2id)
