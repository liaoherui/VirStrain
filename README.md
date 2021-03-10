# VirStrain
An RNA virus strain-level identification tool for short reads.
### E-mail: heruiliao2-c@my.cityu.edu.hk
### Version: V1.0
---------------------------------------------------------------------------
### Dependencies:
* Python >=3.6 (3.6 is recommanded and 3.9 is not supprted now!)
* Perl
* Required python package: networkx, numpy, pandas, biopython, Plotly==3.10.0

(If you have installed conda, then you can run `sh install_package.sh` to install all required packages automatically.)

Make sure these programs have been installed before using VirStrain.

## Install (Linux or ubuntu only)

####
`git clone https://github.com/liaoherui/VirStrain.git`<BR/>
`cd VirStrain`<BR/>
`chmod 755 bin/jellyfish-linux`<BR/>
####

Then, you can download the reference database of 3 RNA viruses. 
There are three ways for you to download the reference database.<BR/><BR/>
-> Method-1:<BR/>
Run:<BR/>
`cd VirStrain`<BR/>
`sh download.sh`<BR/>
-> Method-2:<BR/>
Run:<BR/>
`cd VirStrain`<BR/>
`wget https://github.com/liaoherui/VirStrain/raw/main/VirStrain_DB.tar.gz`<BR/>
`tar -zxvf VirStrain_DB.tar.gz`<BR/>
`rm VirStrain_DB.tar.gz` <BR/>
-> Method-3:<BR/>
If you have installed git lfs, then you can simply run: <BR/>
`git lfs clone https://github.com/liaoherui/VirStrain.git`<BR/>
Then the database will be cloned with the repository.<BR/><BR/>
If all failed, please email to the author to get the database.

## Usage

### Use VirStrain to identify RNA virus strains in short reads.

For SE reads:<BR/>
  `python VirStrain.py -i Test_Data/MT451123_1.fq -d VirStrain_DB/SCOV2 -o MT451123_SE_Test`<BR/>

For PE reads:<BR/>
  `python VirStrain.py -i Test_Data/MT451123_1.fq -p Test_Data/MT451123_2.fq -d VirStrain_DB/SCOV2 -o MT451123_PE_Test`<BR/>

When the virus has high mutation rate, like HIV, you need to add `-m` parameter.

For HIV:<BR/>
  SE reads: `python VirStrain.py -i <Read1> -d VirStrain_DB/HIV -o <Output_dir> -m`<BR/>
  PE reads: `python VirStrain.py -i <Read1> -p <Read2> -d VirStrain_DB/HIV -o <Output_dir> -m`<BR/>

### Use VirStrain to build your own custom database.<BR/>
  `python VirStrain_build.py -i <Input_MSA> -d <Database_Dir>`<BR/>
  
  Note: The format of input MSA should be same as the format of MSA generated by Mafft (https://mafft.cbrc.jp/alignment/software/).<BR/>

## Output Format

The output of VirStrain contains two files. The first is a report file in text format. This file contains all identified strains and their depth and site coverage, etc. The other file is an interactive HTML page to display the depth and uniqueness of sites. 

You can check the output file in the folder "MT451123_Sim_PE".

The picture below displays an output example of a simulated data (Truth: MT451123.1). <BR/>

![VirStrain Report](https://github.com/liaoherui/VirStrain/blob/main/Output_fmt/report_simulate.png)

Explaination about the columns in the output of VirStrain:

Column_name    |	Description	
------------ | ------------- 
Strain_ID |	The NCBI accession number of identified strain.
Cls_info | The cluster information of identified strain, eg: Cluster2830_2 -> belong Cluster2830, size=2.
SubCls_info | The sub-cluster information of identified strain.
Vscore | The Vscore generated by VirStrain algorithm.
Valid_Map_Rate | The covered sites out of total sites in the first iteration of VirStrain.
Site_coverage | The covered sites out of total sites in the current iteration of VirStrain.
Strain_depth | The sequencing depth of identified strain predicted by VirStrain.
Strain_info | The metadata of identified strains, such as region information and subtype, etc.
SNV_freq | The SNV frequency of all sites.



