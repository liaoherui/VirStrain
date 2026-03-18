import re
import os
#import json
import sys
import getopt

input_fa=''
cls_file=''
out_dir=''
opts,args=getopt.getopt(sys.argv[1:],"i:c:o:")
for opt,arg in opts:
	if opt=='-i':
		input_fa=arg
	elif opt=='-c':
		cls_file=arg
	elif opt=='-o':
		out_dir=arg
if not out_dir:
	out_dit='VirStrain_DB'
#f=open('head_file.txt','r')
f=open(input_fa)
dhead={}
while True:
	line=f.readline().strip()
	if not line:break
	if not re.search('>',line):continue
	if re.search('\|',line):
		pre=re.split('\|',line)[0].strip()
	else:
		pre=line.split()[0].strip()
	#print(pre)
	#exit()
	info=re.sub('.*human','',line)
	info=re.sub(',.*','',info)
	dhead[pre]=info

#clade=json.loads('nextclade.json')
'''
with open('nextclade.json','r') as j:
	clade=json.loads(j.read())
dclade={}
for c in clade:
	if 'clade' not in c:
		dclade['>'+c['seqName']]='NA'
	else:
		dclade['>'+c['seqName']]=c['clade']
'''

f2=open(cls_file,'r')
o=open(out_dir+'/Anno_rebuild_cls.clstr','w+')
while True:
	line=f2.readline().strip()
	if not line:break
	if re.search('Cluster',line) or re.search('Cls',line):
		if re.search('Cluster',line):
			o.write(line+'\n')
		else:
			o.write('\t'+line+'\n')
	else:
		ele=line.split('\t')
		o.write('\t\t'+line+'\t'+dhead[ele[0]]+'\n')
