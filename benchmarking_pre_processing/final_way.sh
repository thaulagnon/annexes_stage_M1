#!/bin/bash

# This script is the result of the benchmarking comparing old_way.sh and new_way.sh
# It picks the tool best suited for the job at each step
# And adds a complementary trimming step for unmerged reads, which AdapterRemoval does not trim properly
# It uses AdapterRemoval for trimming and merging
#         bbduk for additional unmerged read trimming
#         Prinseq for complexity filtering
#         bbduk for PhiX filtering
#         Prinseq for duplicate removal

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
echo  `date`

AdapterRemoval --file1 $INPUT_FASTQ_R1 \
               --file2 $INPUT_FASTQ_R2 \
               --settings $OUTPUT_FOLDER/AdapterRemoval_report.txt \
               --output1 $OUTPUT_FOLDER/final_way.trimmed.R1.fq \
               --output2 $OUTPUT_FOLDER/final_way.trimmed.R2.fq \
               --discarded $OUTPUT_FOLDER/AdapterRemoval.discarded.fq \
               --outputcollapsed $OUTPUT_FOLDER/final_way.trimmed.M.fq \
               --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
               --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
               --minlength 25 \
               --collapse \
               --minalignmentlength 15 \
               --qualitymax 45

echo  "Done Trimming and Merging reads."
echo  `date`

##################################################################################################
################################# bbduk : Trim paired end reads ##################################
##################################################################################################

# This is a second trimming step, which applies only to unmerged reads
# AdapterRemoval performs very poorly when trimming these unmerged reads
# Hence we need a second trimming step, using bbduk, whose k-mer approach is much more thorough

BBDUK_ADAPTERS=/home/taulagnon/miniconda3/envs/env_test/share/bbmap/resources/adapters.fa

echo  "Second trimming on" ${OUTPUT_FOLDER}"/new_way.trimmed.R1.fq and" ${OUTPUT_FOLDER}"/new_way.trimmed.R2.fq"
echo  `date`

bbduk.sh in=$OUTPUT_FOLDER/final_way.trimmed.R1.fq \
         in2=$OUTPUT_FOLDER/final_way.trimmed.R2.fq \
         ref=$BBDUK_ADAPTERS \
         out=$OUTPUT_FOLDER/final_way.retrimmed.R1.fq \
         out2=$OUTPUT_FOLDER/final_way.retrimmed.R2.fq \
         ktrim=r \
         k=23 \
         mink=11 \
         minlength=25

echo  "Done re-trimming Unmerged reads"
echo  `date`

##################################################################################################
################################## PRINSEQ : Complexity Filtering ##################################
##################################################################################################

# In the end, seeing how poorly fastp performed, it was decided to stick with PRINSEQ in spite of its old age
# PRINSEQ was implemented into Galaxy, which made it simple to add for our workflow

### Complexity filtering for merged (M) reads

echo "Filtering low-complexity sequences from" ${OUTPUT_FOLDER}"/final_way.trimmed.M.fq"
echo  `date`

prinseq-lite.pl -fastq ${OUTPUT_FOLDER}/final_way.trimmed.M.fq  \
                -out_good ${OUTPUT_FOLDER}/final_way.complexity_filtered.M  \
                -out_bad null  \
                -lc_method dust  \
                -lc_threshold 7  \
                -line_width 0 \
                2>&1

echo  "Done Filtering low-complexity reads from Merged reads"
echo  `date`

### Complexity filtering for unmerged (R1+R2) reads

echo "Filtering low-complexity sequences from" ${OUTPUT_FOLDER}"/final_way.retrimmed.R1.fq and" ${OUTPUT_FOLDER}"/final_way.retrimmed.R2.fq"
echo  `date`

prinseq-lite.pl -fastq ${OUTPUT_FOLDER}/final_way.retrimmed.R1.fq \
                -fastq2 ${OUTPUT_FOLDER}/final_way.retrimmed.R2.fq \
                -out_good ${OUTPUT_FOLDER}/final_way.complexity_filtered \
                -out_bad null \
                -lc_method dust \
                -lc_threshold 7 \
                -line_width 0 \
                2>&1

echo  "Done Filtering low-complexity reads from Unmerged reads"
echo  `date`

### Renaming Prinseq outputs in order to have somewhat consistent names

for i in 1 2
    do
        mv ${OUTPUT_FOLDER}/final_way.complexity_filtered_${i}.fastq ${OUTPUT_FOLDER}/final_way.complexity_filtered.R${i}.fq
    done

##################################################################################################
############################# bbduk : PhiX decontamination #######################################
##################################################################################################

# PhiX filtering is performed because PhiX sequences are spiked in during sequencing
# This is done in order to reach a balanced nucleotide composition in the first 6 nts (Orlando et al 2021)
# Problem is, PhiX DNA gets sequenced and ends up in our data
# Which is a problem for barley, whose genome is poorly assembled and contains PhiX sequences, leading to false mapping.

PHIX_REF=/home/taulagnon/miniconda3/envs/env_test/share/bbmap/resources/phix174_ill.ref.fa.gz

### PhiX filtering for merged (M) reads

echo "Filtering PhiX sequences from" ${OUTPUT_FOLDER}"/final_way.complexity_filtered.M.fq"
echo  `date`

bbduk.sh in=${OUTPUT_FOLDER}/final_way.complexity_filtered.M.fastq \
         out=${OUTPUT_FOLDER}/final_way.phiX_filtered.M.fq \
         ref=$PHIX_REF \
         k=31 hdist=1 qin=33 \
         stats=${OUTPUT_FOLDER}/final_way.phiX_removal_stats.M.txt \
         overwrite=t

echo  "Done filtering PhiX sequences from Merged reads."
echo  `date`

### PhiX filtering for unmerged (R1+R2) reads

echo "Filtering PhiX sequences from" ${OUTPUT_FOLDER}"/final_way.complexity_filtered.R1.fq and" ${OUTPUT_FOLDER}"/final_way.complexity_filtered.R2.fq"
echo  `date`

bbduk.sh in=${OUTPUT_FOLDER}/final_way.complexity_filtered.R1.fq \
         in2=${OUTPUT_FOLDER}/final_way.complexity_filtered.R2.fq \
         out=${OUTPUT_FOLDER}/final_way.phiX_filtered.R1.fq \
         out2=${OUTPUT_FOLDER}/final_way.phiX_filtered.R2.fq \
         ref=$PHIX_REF \
         k=31 hdist=1 qin=33 \
         stats=${OUTPUT_FOLDER}/final_way.phiX_removal_stats.PE.txt \
         overwrite=t

echo  "Done filtering PhiX sequences from Unmerged reads."
echo  `date`

##################################################################################################
############################### fastp : Duplicate Removal ########################################
##################################################################################################

# NB : In the end, this step was not kept in the Galaxy workflow
# Because it was done prior to BLASTing only
# And BLASTing has been abandoned on Galaxy because of resource isssues

### Duplicate Removal for merged (M) reads

echo "Removing duplicates from" ${OUTPUT_FOLDER}"/final_way.phiX_filtered.M.fq"
echo  `date`

fastp -A \
      -L \
      -Q \
      -G \
      -i $OUTPUT_FOLDER/final_way.phiX_filtered.M.fq \
      -o $OUTPUT_FOLDER/final_way.duplicates_removed.M.fq \
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
      -i $OUTPUT_FOLDER/final_way.phiX_filtered.R1.fq \
      -I $OUTPUT_FOLDER/final_way.phiX_filtered.R2.fq \
      -o $OUTPUT_FOLDER/final_way.duplicates_removed.R1.fq \
      -O $OUTPUT_FOLDER/final_way.duplicates_removed.R2.fq \
      --dedup

echo  "Done removing duplicates from Unmerged reads"
echo  `date`
