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

src_path=sys.path[0]

input_msa=''

cls_file=''
dash_cutoff=0.0
k=25
opts,args=getopt.getopt(sys.argv[1:],"i:r:")

for opt,arg in opts:
	if opt=='-i':
		input_msa=arg
	if opt=='-r':
		cls_file=arg	

	

#def cls_with_entropy_hier(alignment, centropy ,cbase_freq, sid, sname2seq, sid2name, sname2id):
def cls_with_perl_hier(sstrain2seq,strain_set,src_path):
	selected_columns={}
	#copyentropy=centropy # copyentropy will be updated during the iterative process
	#copyalignment=alignment # copyalignment will not be updated
	cnum=0
	fs=open('Vbuild_tem_cls.aln','w+')
	for s in sstrain2seq:
		fs.write(s+'\n'+sstrain2seq[s]+'\n')
		cnum=len(sstrain2seq[s])
	fs.close()
	os.system('perl '+src_path+'/aln2cluster-overlap-kmer-withd.pl Vbuild_tem_cls.aln '+str(len(strain_set))+' 0 0 '+str(cnum)+' 0 > Vbuild_tem.column')
	#print('perl aln2cluster-overlap-kmer-withd.pl Vbuild_tem_cls.aln '+str(len(strain_set))+' 0 0 '+str(cnum)+' 0 > Vbuild_tem.column')
	#exit()
	fs2=open('Vbuild_tem.column','r')
	while True:
		line=fs2.readline().strip()
		if not line:break
		if not re.search('column',line):continue
		ele=line.split()
		selected_columns[int(ele[1])]=''
	#os.system('rm Vbuild_tem_cls.aln Vbuild_tem.column')
	#exit()
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
				cseq+=sstrain2seq[s][c]
			tem.append(cseq)
			tem_s=str(strain_set[s])+'_cls'
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
				base=sstrain2seq[s][c].upper()
				left=sstrain2seq[s][:c][::-1]
				right=sstrain2seq[s][c+1:]
				pos_base=str(c)+'-'+base
				if c-12<0:
					block_seq=sstrain2seq[s][0:c]+base+sstrain2seq[s][c+1:c+13]
				else:
					'''
					print(c+13,c-12,base)
					print(alignment[s].seq[c-12:c])
					print(alignment[s].seq[c+1:c+13])
					print(len(alignment[s,:]))
					exit()
					'''
					block_seq=sstrain2seq[s][c-12:c]+base+sstrain2seq[s][c+1:c+13]
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
					kmr_sub[kmr]={strain_sub[strain_set[s]]:''}
				else:
					kmr_sub[kmr][strain_sub[strain_set[s]]]=''
				if kmr not in kmr_pos:
					kmr_pos[kmr]={pos_base:''}
				else:
					kmr_pos[kmr][pos_base]=''
		# recluster according to pos-snp kmer
		#for k in kmr_sub:
		
		return  sub_info,kmr_sub,kmr_pos,sort_sub				
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
	sname2seq['>'+name]=str(record.seq)

#sid2name=dict(zip(sid,sname))  # Strain_ID -> Strain_Name
sid2name={}
sname2id=dict(zip(sname,sid))  # Strain_Name -> Strain_ID

## Check how many columns left
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
			sid2name[int(ele[-1])]=rep
			cls_rep[cls]=rep
		else:
			strain_rep[ele[0]]=rep
			strain_info[ele[0]]=cls
			cls_strain[cls][ele[0]]=int(ele[-1])
			sid2name[int(ele[-1])]=ele[0]

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
			o1.write(s+'\t'+strain_rep[s]+'\t'+strain_info[s]+'\t'+strain_info[s]+'\t'+str(cls_strain[c][s])+'\t1\n')
			o2.write('>'+c+'\n')
			o2.write('\t'+s+'\t*\t'+str(cls_strain[c][s])+'\n')
	else:
		strain_set={}
		ssname2seq={}
		for s in  cls_strain[c]:
			strain_set[s]=cls_strain[c][s]
			ssname2seq[s]=sname2seq[s]
		sub_info,kmr_sub,kmr_pos,sort_sub=cls_with_perl_hier(ssname2seq,strain_set,src_path)
		if len(sub_info)==0 or len(sort_sub)==1:
			new_cls_size.append(len(cls_strain[c]))
			o2.write('>'+c+'\n')
			for s in cls_strain[c]:
				o1.write(s+'\t'+strain_rep[s]+'\t'+strain_info[s]+'\t'+strain_info[s]+'\t'+str(cls_strain[c][s])+'\t'+str(len(cls_strain[c]))+'\n')
				if s==cls_rep[c]:
					o2.write('\t'+s+'\t*\t'+str(cls_strain[c][s])+'\n')
				else:
					o2.write('\t'+s+'\t'+str(cls_strain[c][s])+'\n')
		else:
			o2.write('>'+c+'\n')
			for sub in sort_sub:
				o2.write('\t>'+sub+'\n')
				new_cls_size.append(len(sub_info[sub]))
				for ss in sub_info[sub]:
					ss=sid2name[ss]
					o1.write(ss+'\t'+strain_rep[ss]+'\t'+strain_info[ss]+'\t'+sub+'\t'+str(cls_strain[c][ss])+'\t'+str(len(sub_info[sub]))+'\n')
					if ss==cls_rep[c]:
						o2.write('\t\t'+ss+'\t*\t'+str(cls_strain[c][ss])+'\n')
					else:
						o2.write('\t\t'+ss+'\t'+str(cls_strain[c][ss])+'\n')
			for kl in kmr_sub:
				info1=','.join(list(kmr_pos[kl].keys()))
				info2=','.join(list(kmr_sub[kl].keys()))
				o3.write(kl+'\t'+c+'\t'+info1+'\t'+info2+'\n')
				seq=Seq(kl,IUPAC.ambiguous_dna)
				k2=seq.reverse_complement()._data
				o3.write(k2+'\t'+c+'\t'+info1+'\t'+info2+'\n')
				
os.system('rm Vbuild_tem.column Vbuild_tem_cls.aln')
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
