import re
import os
import getopt
import sys
import networkx
from networkx.algorithms.components.connected import connected_components

input_matrix=''

opts,args=getopt.getopt(sys.argv[1:],"i:")
for opt,arg in opts:
	if opt=='-i':
		input_matrix=arg

f=open(input_matrix,'r')
o1=open('Strain_pos_snp_matrix_not_redundant_MM_Call.txt','w+')
o2=open('Remove_redundant_matrix_MM_Call.clstr','w+')
line=f.readline().strip()
o1.write('\t\t'+line+'\n')
head_line=line.split('\t')
carr=[]
base=['A','T','G','C']
for e in head_line:
	c=re.split('-',e)[0]
	if int(c) not in carr:
		carr.append(int(c))
#o1.write('\t\t'+line+'\n')

arr=[]
ds_snp={}
ds2id={}
#strain=[]
i=0
while True:
	line=f.readline().strip()
	if not line:break
	ele=line.split('\t')
	tem=[]
	tem.append(ele[0])
	ds2id[ele[0]]=i
	i+=1
	#strain.append(ele[0])
	ds_snp[ele[0]]=ele[1:]
	snp_info=','.join(ele[1:])
	tem.append(snp_info)
	arr.append(tem)

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

G=to_graph(arr)
i=1
final_strain=[]
dcsb={} # column -> strain -> base
for c in carr:
	dcsb[c]={}
for ele in list(connected_components(G)):
	ele=list(ele)
	final=[]
	for e in ele:
		if not re.search(',',e):
			final.append(e)
	o2.write('>Cluster'+str(i)+'_'+str(len(final))+'\n')
	i+=1
	o2.write('\t'+final[0]+'\t*\t'+str(ds2id[final[0]])+'\n')
	o1.write(final[0]+'\t'+'\t'.join(ds_snp[final[0]])+'\n')
	final_strain.append(final[0])
	ti=0
	for e in ds_snp[final[0]]:
		if int(e)==1:
			info=head_line[ti]
			arr=re.split('-',info)
			dcsb[int(arr[0])][final[0]]=arr[1]
		ti+=1
	for s in final[1:]:
		o2.write('\t'+s+'\t'+str(ds2id[s])+'\n')
'''
o1.write('POS,')
s1=','.join(final_strain)
o1.write(s1+'\n')
print(len(final_strain))
exit()
for c in carr:
	o1.write(str(c)+',')
	tem=[]
	for s in final_strain:
		tem.append(dcsb[c][s])
	s2=','.join(tem)
	o1.write(s2+'\n')
'''
