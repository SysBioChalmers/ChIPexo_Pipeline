# ChIPexo Pipeline

This repostiory contains the full pipeline to analyse ChIPexo data, starting from raw sequence read files obtained after sequencing.

- Last update: 2019-06-10

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

### Start
The file "ChIPexo_Pipeline_MAIN.bash" is the main script of the pipeline, which manages all calls to other scripts and software packages. It can be called from most terminals which have bash.
In the first part of the script, the options including the paths to the Input data, Temporary storage space and the software dependencies have to be set.
Below one can choose which parts of the pipeline to run, e.g. which output files will be generated.
It is therefore also possible to skip the bowtie2 sequence alignment part if one already has mapped .bam files available. 

## Contributors
- [Christoph Boerlin](https://www.chalmers.se/en/staff/Pages/borlinc.aspx); Chalmers University of Technology, Gothenburg Sweden
- [David Bergenholm](https://www.chalmers.se/en/staff/Pages/david-jullesson.aspx); Chalmers University of Technology, Gothenburg Sweden
- Petter Holland; Chalmers University of Technology, Gothenburg Sweden
