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

Then, you can download the reference database of 3 RNA viruses. Run:<BR/>
`cd VirStrain`<BR/>
`sh download.sh`<BR/>

If failed, please email to the author to get the database.

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

### Use VirStrain to build your custom database.<BR/>
  `python VirStrain_build.py -i <Input_MSA> -d <Database_Dir>`<BR/>
  
  Note: The format of input MSA should be same as the format of MSA generated by Mafft (https://mafft.cbrc.jp/alignment/software/).<BR/>

## Output Format

You can check the output file in the folder "MT451123_Sim_PE".

The picture below displays an output example of a simulated data (Truth: MT451123.1). <BR/>

![VirStrain Report](https://github.com/liaoherui/VirStrain/blob/main/Output_fmt/report_simulate.png)



