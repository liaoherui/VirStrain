import re
import os,sys
import argparse
from bin import prescan
import uuid

__author__="Liao Herui"

usage="VirStrain - An RNA virus strain-level identification tool for short reads."

def split_contig(ingenome,kmers,o):
    # Apply Bowtie2 to index genomes
    uid1=uuid.uuid1().hex
    os.system('bowtie2-build '+ingenome+' '+uid1)
    os.system('bowtie2 -x '+uid1+' -f '+kmers+' --score-min C,0,0 -S '+uid1+'.sam')
    f=open(uid1+'.sam','r')
    dmatch={}
    while True:
        line=f.readline().strip()
        if not line:break
        if line[0]=='@':continue
        ele=line.split('\t')
        if not ele[2]=='\*':
            dmatch[ele[2]]=''

    os.system('rm '+uid1+'.*')
    uid=uuid.uuid1().hex
    if not os.path.exists(uid):
        os.makedirs(uid)
    f=open(ingenome,'r')
    lfilter=100
    d={}
    while True:
        line=f.readline().strip()
        if not line:break
        if re.search('>',line):
            pre=line.split()[0]
            pre=re.sub('>','',pre)
            if pre in dmatch:
                keep=1

            else:
                keep=0
                o.write(pre+'\tNA\tSkip\n')
            tem=line
            if keep==1:
                d[tem]=''
        else:
            if keep==1:
                d[tem]+=line
    res={}
    for s in d:
        if len(d[s])<lfilter:continue
        pre=re.sub('>','',s)
        pre=pre.split()[0]
        o=open(uid+'/'+pre+'.fasta','w+')
        o.write(s+'\n'+d[s]+'\n')
        res[uid+'/'+pre+'.fasta']=''
    #print(res,uid)
    #exit()
    return res,uid

     
def main():
	parser=argparse.ArgumentParser(prog='VirStrain.py',description=usage)
	parser.add_argument('-i','--input_reads',dest='input_reads',type=str,required=True,help="Input fastq data --- Required")

	parser.add_argument('-p','--input_reads2',dest='input_reads2',type=str,help="Input fastq data for PE reads.")
	parser.add_argument('-d','--database_dir',dest='db_dir',type=str,required=True,help='Database dir --- Required')
	parser.add_argument('-o','--output_dir',dest='out_dir',type=str,help='Output dir (default: current dir/VirStrain_Out)')
	parser.add_argument('-c','--site_filter_cutoff',dest='sf_cutoff',type=str,help='The cutoff of filtering one site (default: 0.05)')
	parser.add_argument('-v', '--comprehensive_mode', dest='cv_mode', type=str,help='If set to 1, then VirStrain will identify viral strains for input contigs in a more comprehensive database. (default: 0)')
	parser.add_argument('-s','--rank_by_sites',dest='rk_site',type=str,help='If set to 1, then VirStrain will sort the most possible strain by matches to the sites. (default: 0)')
	parser.add_argument('-m','--high_mutation_virus',dest='hm_virus',help='If the virus has high mutation rate, use this option. (default: Not use)',default='1',nargs='?')

	
	args=parser.parse_args()	
	in_read1=args.input_reads
	in_read2=args.input_reads2
	db_dir=args.db_dir
	out_dir=args.out_dir
	sfc=args.sf_cutoff
	cvm = args.cv_mode
	rks=args.rk_site
	hm=args.hm_virus
	if not cvm:
		cvm=0
	else:
		cvm=int(cvm)
	if not sfc:
		sfc=0.05
	if not rks:
		rks=0
	if not out_dir:
		out_dir='VirStrain_Out'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	if cvm == 1:
		o = open(out_dir + '/VirStrain_contig_report.txt', 'w+')
		o.write('Contigs_ID\tSpecies_info\tStrain_info\n')
		contigs,uid=split_contig(in_read1,db_dir+'/merge_db.fa',o)
		if not os.path.exists(out_dir+'/Report'):
			os.makedirs(out_dir+'/Report')
		for genome in contigs:
			pre=re.split('/',genome)[-1]
			pre=re.sub('\.fasta','',pre)
			o.write(pre)
			target_sp,kmatch=prescan.scan(genome,db_dir+'/merge_db.fa')
			if int(kmatch)<10:
				o.write('\t'+target_sp+':'+str(kmatch)+'\tSkip\n')
				continue
			nd_dir=db_dir+'/'+target_sp
			os.system('python bin/S3_Strain_pred_My_Method_V0819_Val_contig.py -i '+genome+' -s '+nd_dir+'/Pos-snp-kmer-all.txt -m '+nd_dir+'/Strain_pos_snp_matrix_not_redundant_MM_Call.txt -f '+nd_dir+'/Pos-snp-kmer-all.fa -c '+nd_dir+'/Strain_cls_info.txt -b '+nd_dir+'/SubCls_kmer.txt -d '+nd_dir+' -o '+pre+'_VirStrain_report.txt')
			#exit()
			f=open(pre+'_VirStrain_report.txt','r')
			line=f.readline().strip()
			line=f.readline().strip()
			line=f.readline().strip()
			line=re.sub('\t','|',line)
			o.write('\t'+target_sp+':'+str(kmatch)+'\t'+line+'\n')
			os.system('mv '+pre+'_VirStrain_report.txt '+out_dir+'/Report')
			#exit()
		os.system('rm -rf '+uid)
	exit()
	# whether this is the virus with high mutation rate, like HIV
	if not hm:
		if not in_read2:
			os.system('python bin/S3_Strain_pred_My_Method_For_HIV.py -i '+in_read1+' -s '+db_dir+'/Pos-snp-kmer-all.txt  -f '+db_dir+'/Pos-snp-kmer-all.fa -c '+db_dir+'/Remove_redundant_matrix_MM_Call.clstr -d '+db_dir+'/ID2Name.txt  -b '+db_dir+' -o VirStrain_report.txt')
		else:
			os.system('python bin/S3_Strain_pred_My_Method_For_HIV.py -i '+in_read1+' -p '+in_read2+' -s '+db_dir+'/Pos-snp-kmer-all.txt  -f '+db_dir+'/Pos-snp-kmer-all.fa -c '+db_dir+'/Remove_redundant_matrix_MM_Call.clstr -d '+db_dir+'/ID2Name.txt  -b '+db_dir+' -o VirStrain_report.txt')
		if os.path.exists('VirStrain_report.txt'):
			os.system('python bin/S4_Plot_strain_cov.py')
			os.system('mv VirStrain_report.txt VirStrain_report.html Mps_ps_depth.csv Ops_ps_depth.csv '+out_dir)
		exit()
	# Start to identify strains
	if not in_read2: # SE reads
		if rks==0:
			os.system('python bin/S3_Strain_pred_My_Method_V0819_Val.py -i '+in_read1+' -s '+db_dir+'/Pos-snp-kmer-all.txt -m '+db_dir+'/Strain_pos_snp_matrix_not_redundant_MM_Call.txt -f '+db_dir+'/Pos-snp-kmer-all.fa -c '+db_dir+'/Strain_cls_info.txt -b '+db_dir+'/SubCls_kmer.txt -d '+db_dir+' -o VirStrain_report.txt')
		else:
			os.system('python bin/S3_Strain_pred_My_Method_V0819_Val.py -i '+in_read1+' -s '+db_dir+'/Pos-snp-kmer-all.txt -m '+db_dir+'/Strain_pos_snp_matrix_not_redundant_MM_Call.txt -f '+db_dir+'/Pos-snp-kmer-all.fa -c '+db_dir+'/Strain_cls_info.txt -b '+db_dir+'/SubCls_kmer.txt -d '+db_dir+' -r 1 -o VirStrain_report.txt')
	else: # PE reads
		if rks==0:
			os.system('python bin/S3_Strain_pred_My_Method_V0819_Val.py -i '+in_read1+' -p '+in_read2+' -s '+db_dir+'/Pos-snp-kmer-all.txt -m '+db_dir+'/Strain_pos_snp_matrix_not_redundant_MM_Call.txt -f '+db_dir+'/Pos-snp-kmer-all.fa -c '+db_dir+'/Strain_cls_info.txt -b '+db_dir+'/SubCls_kmer.txt -d '+db_dir+' -o VirStrain_report.txt')
		else:
			os.system('python bin/S3_Strain_pred_My_Method_V0819_Val.py -i '+in_read1+' -s '+db_dir+'/Pos-snp-kmer-all.txt -m '+db_dir+'/Strain_pos_snp_matrix_not_redundant_MM_Call.txt -f '+db_dir+'/Pos-snp-kmer-all.fa -c '+db_dir+'/Strain_cls_info.txt -b '+db_dir+'/SubCls_kmer.txt -d '+db_dir+' -r 1 -o VirStrain_report.txt')
	# Plot Html Page
	if os.path.exists('VirStrain_report.txt'):
		os.system('python bin/S4_Plot_strain_cov.py')	
		os.system('mv VirStrain_report.txt VirStrain_report.html Mps_ps_depth.csv Ops_ps_depth.csv '+out_dir)


if __name__=='__main__':
	sys.exit(main())
