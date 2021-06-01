import re
import os,sys
import argparse
import math

__author__="Liao Herui"

usage="VirStrain - An RNA virus strain-level identification tool for short reads."

def main():
	parser=argparse.ArgumentParser(prog='VirStrain_build.py',description=usage)
	parser.add_argument('-i','--input_msa',dest='input_msa',type=str,required=True,help="Input MSA file generated by mafft --- Required")
	parser.add_argument('-d','--database_dir',dest='db_dir',type=str,help='Database Output dir (default: current workdir/VirStrain_DB)')
	parser.add_argument('-c','--dash_cutoff',dest='dash_cutoff',type=str,help='The cutoff of dash in each column of MSA (default: 0)')
	parser.add_argument('-s','--sites_cutoff',dest='sites_cutoff',type=str,help='The cutoff of sites number for manual-covering function ( eg. 1 means all useful sites will be used and 0.8 means 80%% useful sites will be used.')
	
	args=parser.parse_args()	
	in_msa=args.input_msa
	db_dir=args.db_dir
	dashc=args.dash_cutoff
	sitec=args.sites_cutoff
	
	if not db_dir:
		db_dir='VirStrain_DB'
	if not dashc:
		dashc=0
	if not sitec:
		sitec=0
	else:
		sitec=float(sitec)
	if not os.path.exists(db_dir):
		os.makedirs(db_dir)
	# Start to build the database
	snum,cnum=scan_msa(in_msa)
	# Run all program
	# Step1 - Choose sites
	os.system('perl bin/aln2cluster-overlap-kmer-withd.pl '+in_msa+' '+str(snum)+' '+str(dashc)+' 0 '+str(cnum)+' 0 > VirStrain_build.column')
	if sitec>0:
		manual_covering(sitec,snum,in_msa)
	# Step2 - Extract kmer and generate snp matrix
	os.system('python bin/S1_extract_kmer.py -i '+in_msa+' -c VirStrain_build.column')
	# Step3 - Divide strains into clusters
	os.system('python bin/S2_remove_redundant_For_Me.py -i Strain_pos_snp_matrix_consider_all.txt')
	# Step4 - Hierarchical Clustering
	os.system('python bin/S2.5_Split_cls_hierarchical_with_perl.py -i '+in_msa+' -r Remove_redundant_matrix_MM_Call.clstr')
	os.system('mv After_split.png Before_split.png ID2Name.txt '+db_dir)
	os.system('mv Pos-snp-kmer-all.fa Pos-snp-kmer-all.txt Rebuild_cls.clstr Remove_redundant_matrix_MM_Call.clstr  Strain_pos_snp_matrix_not_redundant_MM_Call.txt Strain_pos_snp_matrix_consider_all.txt SubCls_kmer.txt Strain_cls_info.txt VirStrain_build.column '+db_dir)
	
def manual_covering(sitec,snum,in_msa):
	os.system('perl bin/aln2entropy.pl '+in_msa+' '+str(snum)+' > manual.column')
	#print('perl bin/aln2entropy.pl '+in_msa+' '+str(snum)+'> manual.column')
	#exit()
	f=open('manual.column','r')
	o=open('mc.column','w+')
	a=[]
	while True:
		line=f.readline().strip()
		if not line:break
		if not re.search('column',line):continue
		ele=line.split()
		for e in ele:
			if re.search('-',e):
				dash_num=re.sub('-:','',e)
				break
		dash_num=int(dash_num)
		if dash_num==0:
			a.append(line)
			#o.write(line+'\n')
	out_c=math.ceil(sitec*float(len(a)))
	#out_c=int(out_c)
	outa=a[:out_c]
	print('Newly added sites number: ',len(outa),'/',len(a))
	for line in outa:
		o.write(line+'\n')
	o.close()
	os.system('rm manual.column')
	os.system('cat VirStrain_build.column mc.column > VirStrain_build_tem.column')
	os.system('rm VirStrain_build.column mc.column')
	os.system('mv VirStrain_build_tem.column VirStrain_build.column')

def scan_msa(in_msa):
	fs1=open(in_msa,'r')
	dseq={}
	while True:
		line=fs1.readline().strip()
		if not line:break
		if re.search('>',line):
			name=line
			dseq[name]=''
		else:
			dseq[name]+=line
	cnum=0
	for d in dseq:
		cnum=len(dseq[d])
		break
	return len(dseq),cnum

if __name__=='__main__':
	sys.exit(main())
