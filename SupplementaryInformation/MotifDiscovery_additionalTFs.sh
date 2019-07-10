#!/bin/bash
#
#Run motif discovery for TFs shown in the SI
#@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
#

mainPath=$(pwd)
outputPath=${mainPath}/SupplementaryInformation
pythonPath=${mainPath}/PythonScripts
refGenomeName=CENPK113-7D
refGenomeBowtiePath=${mainPath}/Data/RefGenomeBowtie

TF=Stb5
cond=Nit
echo "$(date +%T) ${TF} ${cond} extract peak sequences"
python3 ${pythonPath}/extractPeakSequences.py ${TF} ${cond} "${outputPath}/GEM_Files/${TF}_GEM.GEM_events.txt" "${outputPath}/${TF}_${cond}_PeakSequences.bed"
bedtools getfasta -fi ${refGenomeBowtiePath}/${refGenomeName}.fasta -bed ${outputPath}/${TF}_${cond}_PeakSequences.bed -fo ${outputPath}/${TF}_${cond}_PeakSequences.fasta
echo "$(date +%T) ${TF}_${cond} run MEME"
meme ${outputPath}/${TF}_${cond}_PeakSequences.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 5 -maxw 20 -objfun classic -revcomp -markov_order 0 -o ${outputPath}/${TF}_${cond}_MEME

TF=Gcn4
cond=Glu
echo "$(date +%T) ${TF} ${cond} extract peak sequences"
python3 ${pythonPath}/extractPeakSequences.py ${TF} ${cond} "${outputPath}/GEM_Files/${TF}_GEM.GEM_events.txt" "${outputPath}/${TF}_${cond}_PeakSequences.bed"
bedtools getfasta -fi ${refGenomeBowtiePath}/${refGenomeName}.fasta -bed ${outputPath}/${TF}_${cond}_PeakSequences.bed -fo ${outputPath}/${TF}_${cond}_PeakSequences.fasta
echo "$(date +%T) ${TF}_${cond} run MEME"
meme ${outputPath}/${TF}_${cond}_PeakSequences.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 5 -maxw 20 -objfun classic -revcomp -markov_order 0 -o ${outputPath}/${TF}_${cond}_MEME

TF=Cbf1
cond=Ana
echo "$(date +%T) ${TF} ${cond} extract peak sequences"
python3 ${pythonPath}/extractPeakSequences.py ${TF} ${cond} "${outputPath}/GEM_Files/${TF}_GEM.GEM_events.txt" "${outputPath}/${TF}_${cond}_PeakSequences.bed"
bedtools getfasta -fi ${refGenomeBowtiePath}/${refGenomeName}.fasta -bed ${outputPath}/${TF}_${cond}_PeakSequences.bed -fo ${outputPath}/${TF}_${cond}_PeakSequences.fasta
echo "$(date +%T) ${TF}_${cond} run MEME"
meme ${outputPath}/${TF}_${cond}_PeakSequences.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 5 -maxw 20 -objfun classic -revcomp -markov_order 0 -o ${outputPath}/${TF}_${cond}_MEME

