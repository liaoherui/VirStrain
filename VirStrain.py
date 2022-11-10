import re
import os,sys
import argparse

__author__="Liao Herui"

usage="VirStrain - An RNA virus strain-level identification tool for short reads."

def main():
	parser=argparse.ArgumentParser(prog='VirStrain.py',description=usage)
	parser.add_argument('-i','--input_reads',dest='input_reads',type=str,required=True,help="Input fastq data --- Required")
	parser.add_argument('-p','--input_reads2',dest='input_reads2',type=str,help="Input fastq data for PE reads.")
	parser.add_argument('-d','--database_dir',dest='db_dir',type=str,required=True,help='Database dir --- Required')
	parser.add_argument('-o','--output_dir',dest='out_dir',type=str,help='Output dir (default: current dir/VirStrain_Out)')
	parser.add_argument('-c','--site_filter_cutoff',dest='sf_cutoff',type=str,help='The cutoff of filtering one site (default: 0.05)')
	parser.add_argument('-s','--rank_by_sites',dest='rk_site',type=str,help='If set to 1, then VirStrain will sort the most possible strain by matches to the sites. (default: 0)')
	parser.add_argument('-m','--high_mutation_virus',dest='hm_virus',help='If the virus has high mutation rate, use this option. (default: Not use)',default='1',nargs='?')

	
	args=parser.parse_args()	
	in_read1=args.input_reads
	in_read2=args.input_reads2
	db_dir=args.db_dir
	out_dir=args.out_dir
	sfc=args.sf_cutoff
	rks=args.rk_site
	hm=args.hm_virus
	if not sfc:
		sfc=0.05
	if not rks:
		rks=0
	if not out_dir:
		out_dir='VirStrain_Out'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
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
