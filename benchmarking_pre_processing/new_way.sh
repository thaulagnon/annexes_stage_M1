#!/bin/bash

# This script pre-processes data with some new tools, with the objective of modernizing the pipeline
# These tools have been chosen because of their ease of access on Galaxy, their modernity, and their widespread use
# It uses AdapterRemoval for trimming and merging
#         fastp for complexity filtering
#         bbduk for PhiX filtering
#         fastp for duplicate removal

##################################################################################################
########################################## Script usage ##########################################
##################################################################################################

usage() { echo "Usage: $0 -1 <fastq file> -2 <fastq file> -f <ouput folder>" 1>&2; exit 1; }

while getopts ":1:2:f:" opt; do
    case $opt in
        1) INPUT_FASTQ_R1=$OPTARG
        ;;
        2) INPUT_FASTQ_R2=$OPTARG
        ;;
        f) OUTPUT_FOLDER=$OPTARG
    esac
done

if  [ -z $INPUT_FASTQ_R1 ] || [ -z $INPUT_FASTQ_R2 ] || [ -z $OUTPUT_FOLDER ]
then
 usage
fi

mkdir -p $OUTPUT_FOLDER

##################################################################################################
############################### AdapterRemoval : Trimming; Merging ###############################
##################################################################################################

# May add quality-based trimming later on, so far, just trying to get an equivalent to SeqPrep2 trimming so as to compare methods

echo  "Trimming and Merging" ${INPUT_FASTQ_R1} "and" ${INPUT_FASTQ_R2}
date=`date`
echo  ${date}

AdapterRemoval --file1 $INPUT_FASTQ_R1 \
               --file2 $INPUT_FASTQ_R2 \
               --settings $OUTPUT_FOLDER/AdapterRemoval_report.txt \
               --output1 $OUTPUT_FOLDER/new_way.trimmed.R1.fq \
               --output2 $OUTPUT_FOLDER/new_way.trimmed.R2.fq \
               --outputcollapsed $OUTPUT_FOLDER/new_way.trimmed.M.fq \
               --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
               --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
               --minlength 25 \
               --collapse \
               --minalignmentlength 15 \
               --qualitymax 45

echo  "Done Trimming and Merging reads."
echo  `date`

##################################################################################################
################################## fastp : Complexity Filtering ##################################
##################################################################################################

# Not happy with fastp's complexity functionalities.
# Defines complexity in a simplistic way (% of bases identical to next bases) which misses a lot of low-complexity sequences
# e.g. GTGTGTGT is a perfectly fine, complex sequence according to fastp
# Which also means it's hard to make it equivalent to Prinseq in the old way in terms of what it does
# Which is why I've left the complexity threshold at its default value

### Complexity filtering for merged (M) reads

echo "Filtering low-complexity sequences from" ${OUTPUT_FOLDER}"/new_way.trimmed.M.fq"
echo  `date`

fastp -A \
      -L \
      -Q \
      -G \
      -i $OUTPUT_FOLDER/new_way.trimmed.M.fq \
      -o $OUTPUT_FOLDER/new_way.complexity_filtered.M.fq \
      -y

echo  "Done Filtering low-complexity reads from Merged reads"
echo  `date`

### Complexity filtering for unmerged (R1+R2) reads

echo "Filtering low-complexity sequences from" ${OUTPUT_FOLDER}"/new_way.trimmed.R1.fq and" ${OUTPUT_FOLDER}"/new_way.trimmed.R2.fq"

fastp -A \
      -L \
      -Q \
      -G \
      -i $OUTPUT_FOLDER/new_way.trimmed.R1.fq \
      -I $OUTPUT_FOLDER/new_way.trimmed.R2.fq \
      -o $OUTPUT_FOLDER/new_way.complexity_filtered.R1.fq \
      -O $OUTPUT_FOLDER/new_way.complexity_filtered.R2.fq \
      -y

echo  "Done Filtering low-complexity reads from Unmerged reads"
echo  `date`

##################################################################################################
############################# bbduk : PhiX decontamination #######################################
##################################################################################################

PHIX_REF=/home/taulagnon/miniconda3/envs/env_test/share/bbmap/resources/phix174_ill.ref.fa.gz

### PhiX filtering for merged (M) reads

echo "Filtering PhiX sequences from" ${OUTPUT_FOLDER}"/new_way.complexity_filtered.M.fq"
echo  `date`

bbduk.sh in=${OUTPUT_FOLDER}/new_way.complexity_filtered.M.fq \
         out=${OUTPUT_FOLDER}/new_way.phiX_filtered.M.fq \
         ref=$PHIX_REF \
         k=31 hdist=1 qin=33 \
         stats=${OUTPUT_FOLDER}/new_way.phiX_removal_stats.M.txt \
         overwrite=t

echo  "Done filtering PhiX sequences from Merged reads."
echo  `date`

### PhiX filtering for unmerged (R1+R2) reads

echo "Filtering PhiX sequences from" ${OUTPUT_FOLDER}"/new_way.complexity_filtered.R1.fq and" ${OUTPUT_FOLDER}"/new_way.complexity_filtered.R2.fq"
echo  `date`

bbduk.sh in=${OUTPUT_FOLDER}/new_way.complexity_filtered.R1.fq \
         in2=${OUTPUT_FOLDER}/new_way.complexity_filtered.R2.fq \
         out=${OUTPUT_FOLDER}/new_way.phiX_filtered.R1.fq \
         out2=${OUTPUT_FOLDER}/new_way.phiX_filtered.R2.fq \
         ref=$PHIX_REF \
         k=31 hdist=1 qin=33 \
         stats=${OUTPUT_FOLDER}/new_way.phiX_removal_stats.PE.txt \
         overwrite=t

echo  "Done filtering PhiX sequences from Unmerged reads."
echo  `date`

##################################################################################################
############################### fastp : Duplicate Removal ########################################
##################################################################################################

### Duplicate Removal for merged (M) reads

echo "Removing duplicates from" ${OUTPUT_FOLDER}"/new_way.phiX_filtered.M.fq"
echo  `date`

fastp -A \
      -L \
      -Q \
      -G \
      -i $OUTPUT_FOLDER/new_way.phiX_filtered.M.fq \
      -o $OUTPUT_FOLDER/new_way.duplicates_removed.M.fq \
      --dedup

echo  "Done removing duplicates from Merged reads"
echo  `date`

### Duplicate Removal for unmerged (R1+R2) reads

echo "Removing duplicates from" ${OUTPUT_FOLDER}"/new_way.phiX_filtered.R1.fq and" ${OUTPUT_FOLDER}"/new_way.phiX_filtered.R2.fq"
echo  `date`

fastp -A \
      -L \
      -Q \
      -G \
      -i $OUTPUT_FOLDER/new_way.phiX_filtered.R1.fq \
      -I $OUTPUT_FOLDER/new_way.phiX_filtered.R2.fq \
      -o $OUTPUT_FOLDER/new_way.duplicates_removed.R1.fq \
      -O $OUTPUT_FOLDER/new_way.duplicates_removed.R2.fq \
      --dedup

echo  "Done removing duplicates from Unmerged reads"
echo  `date`
