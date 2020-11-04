import re
import os
import networkx
from networkx.algorithms.components.connected import connected_components

f=open('H1N1-16273-one-line.fasta','r')

arr=[]
dseq={}
while True:
	line=f.readline().strip()
	if not line:break
	if re.search('>',line):
		name=line
		dseq[name]=''
	else:
		dseq[name]=line

for s in dseq:
	arr.append([s,dseq[s]])

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
o=open('H1N1-16273-one-line-removeSame.fasta','w+')
o2=open('H1N1-16273-one-line-removeSame.clstr','w+')
for ele in list(connected_components(G)):
	ele=list(ele)
	tem=[]
	for e in ele:
		if re.search('>',e):
			tem.append(e)
	o.write(tem[0]+'\n'+dseq[tem[0]]+'\n')
	tem2=[]
	for t in tem:
		l=len(dseq[t])
		s=t+':'+str(l)
		tem2.append(s)
	outs='\t'.join(tem2)
	o2.write(outs+'\n')
