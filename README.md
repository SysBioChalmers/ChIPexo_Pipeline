# ChIPexo Pipeline

This repostiory contains the full pipeline to analyse ChIPexo data, starting from raw sequence read files obtained after sequencing.

- Last update: 2019-07-10

This repository is administered by Christoph Boerlin (https://github.com/ChristophBoerlin), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology

## Using the pipeline

### Required Software
* Python 3.x , including the following packages:
	* numpy
	* scipy
	* pandas
	* matplotlib

* Bowtie2  (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)   Version used 2.3.4.1
* Samtools (http://www.htslib.org/)                                  Version used 1.7
* Bedtools (https://bedtools.readthedocs.io/en/latest/index.html#)   Version used 2.26.0
* BamUtils (https://github.com/statgen/bamUtil)
* GEM      (http://groups.csail.mit.edu/cgs/gem/)                    Version used 3.4
* MEME     (http://meme-suite.org/tools/meme)                        Version used 5.0.5

### How to setup the pipeline to run it on the provided example data:
1) Clone / download the code from Github
2) Download the Ino2 data from Zenodo (https://doi.org/10.5281/zenodo.3242510) in to the Data folder, filename should be Ino2_rawdata.tar.gz.
3) In the script update the paths to Bowtie2, samtools, bedtools and bamUtils so that they can be added to the $PATH linux variable, if they are not already added and callable (lines 26 to 30).
4) Set the paths for the output and data files if needed (lines 13 to 23). The basic setup is using relative path that do not need to be changed.
5) Select which parts of the pipeline should be run (lines 40 to 48) and how many cores are available for computation (line 50).
6) To run it, use a Linux terminal, navigate to the main folder where the ChIPexo_Pipeline_MAIN.bash file is located and execute it using “bash ChIPexo_Pipeline_MAIN.bash”.

### How to run it for a different Transcription factor or condition:
1) Change the variable TF in line 33 and / or variable condList in line 37 of the ChipExo_Pipeline_MAIN.bash file
2) Create a file called TFName_sequenceLength.txt to replace Ino2_sequenceLength.txt in the Data folder specifying the sequence length of the chosen TF. The file should contain one line “seqLength=X” where X is the actual sequence length in bp of the chosen TF.
3) Create a file called TFName_seqFiles.txt to replace Ino2_seqFiles.txt in the Data folder. In this file the names of the sequencing files for each sample should be specified for Read 1 and Read 2. If multiple files should be used one can separate the files using a colon. (For example: seqFiles[‘NewTF_NewCondition_R1’]=”File1_R1.fastq.gz,File2_R1.fastq.gz” )

### How to change to another organism / genome version:
In order to run the pipeline for another organism or genome version the following files are needed to replace the CENPK113-7D files distributed with the pipeline:
1) Update the variable refGenomeName in line 22 of the ChipExo_Pipeline_MAIN.bash file to the name of the new reference genome
2) Genome sequence as a fasta file to replace file CENPK113-7D.fasta in Data/RefGenomeBowtie
3) Bowtie2 index of the genome to replace CENPK113-7D.XXX files in Data/RefGenomeBowtie
4) Sequence of each individual chromosome as a fasta file starting with chr for GEM, replacing the files in Data/RefGenome
5) Text file with the length of each chromosome, replacing CENPK113-7D_chromSizes in Data/RefGenome
6) If GEM should filter out a certain region of the genome, this needs to be specified in the file GEMexclude.txt in Data/RefGenome
7) If the pipeline should filter out a certain region of the genome or whole chromosomes (e.g. the mitochondria) this needs to be specified in the files Filterlist_regions.txt and Filterlist_chromosomes.txt in the Data folder.
8) Transcription Start Site annotations to replace TSSData.tsv in the Data folder. This file should be tab separated and list a transcription start site including the strand for each gene. It has to include these four columns in this order: Name of the gene, Chromosome, TSS position, Strand.


## Contributors
- [Christoph Boerlin](https://www.chalmers.se/en/staff/Pages/borlinc.aspx); Chalmers University of Technology, Gothenburg Sweden
- [David Bergenholm](https://www.chalmers.se/en/staff/Pages/david-jullesson.aspx); Chalmers University of Technology, Gothenburg Sweden
- Petter Holland; Chalmers University of Technology, Gothenburg Sweden
