#!/bin/bash
#
#Main script of the ChIPexo data analysis pipeline
#@author: Christoph S. Boerlin; Chalmers University of Technology, Gothenburg Sweden
#

###################################################################
###################################################################
###################### Options setting part #######################
###################################################################
###################################################################

# Set Paths
mainPath=$(pwd)
outputPath=${mainPath}/ExampleResults
tmpPath=${mainPath}/TMP
softwarePath=${mainPath}/3rdPartySoftware
pythonPath=${mainPath}/PythonScripts
dataPath=${mainPath}/Data
rawDataPath=${dataPath}
refGenomePath=${dataPath}/RefGenome
refGenomeName=CENPK113-7D
refGenomeBowtiePath=${dataPath}/RefGenomeBowtie

#Add bowtie2, samtools, bedtools and bamutils to path if needed
#PATH=${softwarePath}/bowtie2-2.3.3.1/:$PATH
#PATH=${softwarePath}/samtools-1.6/:$PATH
#PATH=${softwarePath}/bedtools2/bin/:$PATH
#PATH=${softwarePath}/Meme/bin/:$PATH
PATH=${softwarePath}/bamUtil-master/:$PATH

#Set TF and Date (used as a postfix for result files)
TF=Ino2
date=190710

#Set names of conditions and replicates
condList=(Eth Glu)
repNames=(1 2)

# Choose which parts of the programm should be run by setting their value to 1 (turn off with 0).
mapFastq=1
bamOut=1
strandSepWigOut=1
overlapWigOut=1
runGEM=1
runAnalysisGEM=1
runAnalysisReads=1
runMotifDiscovery=1

#Set number of available cores, used for bowtie2
numCores=4

###################################################################
###################################################################
####################### Code execution part #######################
###################################################################
###################################################################

# #fix perl language warnings
# export LANGUAGE=en_US.UTF-8
# export LC_ALL=en_US.UTF-8
# export LANG=en_US.UTF-8
# export LC_TYPE=en_US.UTF-8

#create list of all samples
condListWithReps=()
for cond in ${condList[@]}; do
	for repI in ${repNames[@]}; do
		condListWithReps+=(${cond}${repI})
	done
done

#Load sequence length and calculate trim length
source ${rawDataPath}/${TF}_sequenceLength.txt
trim=$(python3 ${pythonPath}/calculateTrimLength.py ${seqLength})

echo ${TF} trim ${trim}

#Load file names for sequence files
source ${rawDataPath}/${TF}_seqFiles.txt

if [ ${mapFastq} == 1 ]; then
	echo "$(date +%T) Unpack tarball to TMPDIR"
	tar xzvf ${mainPath}/Data/${TF}*rawdata.tar.gz -C ${tmpPath}/

	for cond in ${condListWithReps[@]}; do
		echo "$(date +%T) Map reads with bowtie2 for $cond"
		bowtie2 -p ${numCores} -x ${refGenomeBowtiePath}/${refGenomeName} -1 ${tmpPath}/${seqFiles[${cond}_R1]} -2 ${tmpPath}/${seqFiles[${cond}_R2]} -S ${tmpPath}/${TF}_${cond}.sam --no-mixed --no-discordant
	done

	for i in ${condListWithReps[@]}; do
		echo "$(date +%T) ${TF}_${i} generate .bam"
		samtools view -b -u -q 20 -@ $((${numCores}-1)) ${tmpPath}/${TF}_${i}.sam | samtools sort  -n - -o ${tmpPath}/${TF}_${i}_nameSort.bam
		samtools fixmate -m ${tmpPath}/${TF}_${i}_nameSort.bam  ${tmpPath}/${TF}_${i}_fix.bam
		samtools sort ${tmpPath}/${TF}_${i}_fix.bam -o ${tmpPath}/${TF}_${i}_fixSort.bam
		samtools markdup -r ${tmpPath}/${TF}_${i}_fixSort.bam  ${tmpPath}/${TF}_${i}_remDup.bam
		samtools view -b -@ $((${numCores}-1)) -f 0x40 ${tmpPath}/${TF}_${i}_remDup.bam > ${tmpPath}/${TF}_${i}.bam
		samtools index ${tmpPath}/${TF}_${i}.bam
	done

	if [ ${bamOut} == 1 ]; then
		for i in ${condListWithReps[@]}; do
			cp ${tmpPath}/${TF}_${i}.bam* ${outputPath}/
		done
	fi
fi

if [ ${strandSepWigOut} == 1 ]; then
	for i in ${condListWithReps[@]}; do
		echo "$(date +%T) ${TF}_${i} generate stranded read counts trim 1 and combine replicates"
		bedtools genomecov -ibam ${tmpPath}/${TF}_${i}.bam -fs 1 -d -strand + | awk 'NR==1{print "track type=track1"}{if($3>0 && $2!=last+1){print "variableStep chrom="$1"\n"$2,$3;last=$2} else if($3>0 && $2==last+1){print $2,$3;last=$2}}' > ${tmpPath}/${TF}_${i}_singlePos_plus.wig
		bedtools genomecov -ibam ${tmpPath}/${TF}_${i}.bam -fs 1 -d -strand - | awk 'NR==1{print "track type=track1"}{if($3>0 && $2!=last+1){print "variableStep chrom="$1"\n"$2,$3;last=$2} else if($3>0 && $2==last+1){print $2,$3;last=$2}}' > ${tmpPath}/${TF}_${i}_singlePos_minus.wig
	done
	for i in ${condList[@]}; do
		echo "$(date +%T) ${TF}_${i} combine replicates"
		python3 ${pythonPath}/combineReplicateWigFiles.py "${tmpPath}/${TF}_${i}1_singlePos_plus.wig" "${tmpPath}/${TF}_${i}2_singlePos_plus.wig" "${outputPath}/${TF}_${i}_plus_singlePos_combRep_${date}.wig" "${dataPath}/Filterlist_regions.txt" "${dataPath}/Filterlist_chromosomes.txt"
		python3 ${pythonPath}/combineReplicateWigFiles.py "${tmpPath}/${TF}_${i}1_singlePos_minus.wig" "${tmpPath}/${TF}_${i}2_singlePos_minus.wig" "${outputPath}/${TF}_${i}_minus_singlePos_combRep_${date}.wig" "${dataPath}/Filterlist_regions.txt" "${dataPath}/Filterlist_chromosomes.txt"
	done
fi

if [ ${overlapWigOut} == 1 ]; then
	for i in ${condListWithReps[@]}; do
		echo "$(date +%T) ${TF}_${i} generate stranded read counts for trimlength ${trim}"
		bedtools genomecov -ibam ${tmpPath}/${TF}_${i}.bam -fs ${trim} -d -strand + > ${tmpPath}/${TF}_${i}_plus
		bedtools genomecov -ibam ${tmpPath}/${TF}_${i}.bam -fs ${trim} -d -strand - > ${tmpPath}/${TF}_${i}_minus
		echo "$(date +%T) ${TF}_${i} combine overlap"
		paste ${tmpPath}/${TF}_${i}_plus ${tmpPath}/${TF}_${i}_minus -d "\t" | awk -F'\t' 'OFS="\t" {if ($3>$6){print $1,$2,$6}}{if ($3<=$6){print $1,$2,$3}}' | awk 'NR==1{print "track type=track1"}{if($3>0 && $2!=last+1){print "variableStep chrom="$1"\n"$2,$3;last=$2} else if($3>0 && $2==last+1){print $2,$3;last=$2}}' > ${tmpPath}/${TF}_${i}_ol.wig
	done
	for i in ${condList[@]}; do
		echo "$(date +%T) ${TF}_${i} combine replicates"
		python3 ${pythonPath}/combineReplicateWigFiles.py "${tmpPath}/${TF}_${i}1_ol.wig" "${tmpPath}/${TF}_${i}2_ol.wig" "${outputPath}/${TF}_${i}_ol_combRep_${date}.wig" "${dataPath}/Filterlist_regions.txt" "${dataPath}/Filterlist_chromosomes.txt"
		echo "$(date +%T) ${TF}_${i} assign data to genes"
		python3 ${pythonPath}/assignWigDataToGenes.py "${dataPath}/TSSdata.tsv" "${outputPath}/${TF}_${i}_ol_combRep_${date}.wig" "${outputPath}/${TF}_${i}_ol_combRep_geneAssigned_${date}.wigLike"
	done
fi

if [ ${runGEM} == 1 ]; then
	((trimBam = 75 - ${trim}))
	for i in ${condListWithReps[@]}; do
		echo "$(date +%T) ${TF}_${i} creating trimmed .bam ${trim}"
		bam trimBam ${tmpPath}/${TF}_${i}.bam ${tmpPath}/${TF}_${i}_${trim}.bam -R ${trimBam} --clip
		samtools sort -@10 ${tmpPath}/${TF}_${i}_${trim}.bam  -o ${tmpPath}/${TF}_${i}_${trim}.sorted.bam
		samtools index ${tmpPath}/${TF}_${i}_${trim}.sorted.bam
	done
	gemGenomeInput=" --d ${softwarePath}/GEM/Read_Distribution_ChIP-exo.txt --g ${refGenomePath}/${refGenomeName}_chromSizes --genome ${refGenomePath}/ --ex ${refGenomePath}/GEMexclude.txt"
	gemDataInput=''
	for cond in ${condList[@]}; do
		for repI in ${repNames[@]}; do
			gemDataInput+=" --exptCond${cond} ${tmpPath}/${TF}_${cond}${repI}_${trim}.sorted.bam"
		done
	done
	java -jar ${softwarePath}/GEM/gem34.jar ${gemGenomeInput} ${gemDataInput} --f SAM --out ${outputPath}/${TF}_GEM/ --q 2 --k_min 5 --k_max 18 --smooth 3 --min 5 --mrc 50
fi

if [ ${runAnalysisReads} == 1 ]; then
	echo "$(date +%T) ${TF} create PairwiseComparision"
	python3 ${pythonPath}/plotSampleCorrelation.py ${TF} ${tmpPath}/${TF}_SAMPLE_ol.wig ${outputPath} ${refGenomePath}/${refGenomeName}_chromSizes ${date} "${dataPath}/Filterlist_regions.txt" "${dataPath}/Filterlist_chromosomes.txt" ${condListWithReps[@]}
	for cond in ${condList[@]}; do
		echo "$(date +%T) ${TF}_${cond} create readProfile"
		python3 ${pythonPath}/plotTFReadProfile.py ${TF} ${cond} "${outputPath}/${TF}_${cond}_ol_combRep_geneAssigned_${date}.wigLike" "${outputPath}" "${date}"
	done
fi
if [ ${runAnalysisGEM} == 1 ]; then
	echo "$(date +%T) ${TF} analyse GEM results"
	python3 ${pythonPath}/mapGEMpeaks.py ${TF} "${outputPath}/${TF}_GEM/${TF}_GEM.GEM_events.txt" "${dataPath}/TSSdata.tsv" "${outputPath}" "${date}"
	for cond in ${condList[@]}; do
		echo "$(date +%T) ${TF}_${cond} create peakCenteredFigures"
		python3 ${pythonPath}/plotPeakCenteredFigures.py ${TF} ${cond} "${outputPath}/${TF}_GEM/${TF}_GEM.GEM_events.txt" "${outputPath}/${TF}_${cond}_STRAND_singlePos_combRep_${date}.wig" "${outputPath}/${TF}_${cond}" "${date}"
	done
fi
if [ ${runMotifDiscovery} == 1 ]; then
	echo "$(date +%T) ${TF} discover Motifs"
	for cond in ${condList[@]}; do
		echo "$(date +%T) ${TF}_${cond} extract peak sequences"
		python3 ${pythonPath}/extractPeakSequences.py ${TF} ${cond} "${outputPath}/${TF}_GEM/${TF}_GEM.GEM_events.txt" "${outputPath}/${TF}_${cond}_PeakSequences.bed"
		bedtools getfasta -fi ${refGenomeBowtiePath}/${refGenomeName}.fasta -bed ${outputPath}/${TF}_${cond}_PeakSequences.bed -fo ${outputPath}/${TF}_${cond}_PeakSequences.fasta
		echo "$(date +%T) ${TF}_${cond} run MEME"
		meme ${outputPath}/${TF}_${cond}_PeakSequences.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 5 -maxw 20 -objfun classic -revcomp -markov_order 0 -o ${outputPath}/${TF}_${cond}_MEME
	done
fi
echo "$(date +%T) Done"
