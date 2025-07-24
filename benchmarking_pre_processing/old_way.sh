#!/bin/bash

##################################################################################################
########################################## Script usage ##########################################
##################################################################################################

usage() { echo "Usage: $0 -1 <forward fastq file> -2 <reverse fastq file> -f <ouput folder>" 1>&2; exit 1; }

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
################################## SeqPrep2 : Trimming; Merging ##################################
##################################################################################################

ADAPTER_1="AGATCGGAAGAGCACACGTC"
ADAPTER_2="AGATCGGAAGAGCGTCGTGT"

SEQPREP_MIN_LENGTH=25        # Removes unmappably short reads.
SEQPREP_OVERLAP=15           # Allows for confident merging of reads. Can be reduced to 10 if needed.
SEQPREP_QUALITY_CUTOFF=15    # Quality score cutoff for mismatches to be counted in overlap

echo "Trimming and Merging" ${INPUT_FASTQ_R1} "and" ${INPUT_FASTQ_R2}

SeqPrep2 -f $INPUT_FASTQ_R1 \
         -r $INPUT_FASTQ_R2 \
         -1 ${OUTPUT_FOLDER}/old_way.trimmed.R1.fq.gz -2 ${OUTPUT_FOLDER}/old_way.trimmed.R2.fq.gz \
         -q $SEQPREP_QUALITY_CUTOFF \
         -L $SEQPREP_MIN_LENGTH \
         -A $ADAPTER_1 \
         -B $ADAPTER_2 \
         -s ${OUTPUT_FOLDER}/old_way.trimmed.M.fq.gz \
         -E ${OUTPUT_FOLDER}/old_way_readable_alignment.txt.gz \
         -o ${SEQPREP_OVERLAP} \
         -d 1 \
         -C ATCTCGTATGCCGTCTTCTGCTTG \
         -D GATCTCGGTGGTCGCCGTATCATT \
         2>&1

## gunzipping results to allow their treatment by Prinseq

gunzip ${OUTPUT_FOLDER}/old_way.trimmed.R1.fq.gz 
gunzip ${OUTPUT_FOLDER}/old_way.trimmed.R2.fq.gz
gunzip ${OUTPUT_FOLDER}/old_way.trimmed.M.fq.gz

##################################################################################################
############################# Prinseq : Complexity filtering #####################################
##################################################################################################

### Complexity filtering for merged (M) reads

echo "Filtering low-complexity sequences from" ${OUTPUT_FOLDER}"/old_way.trimmed.M.fq"

prinseq-lite.pl -fastq ${OUTPUT_FOLDER}/old_way.trimmed.M.fq  \
                -out_good ${OUTPUT_FOLDER}/old_way.complexity_filtered.M  \
                -out_bad null  \
                -lc_method dust  \
                -lc_threshold 7  \
                -line_width 0 \
                2>&1 

### Complexity filtering for unmerged (R1+R2) reads

echo "Filtering low-complexity sequences from" ${OUTPUT_FOLDER}"/old_way.trimmed.R1.fq and" ${OUTPUT_FOLDER}"/old_way.trimmed.R2.fq"

prinseq-lite.pl -fastq ${OUTPUT_FOLDER}/old_way.trimmed.R1.fq \
                -fastq2 ${OUTPUT_FOLDER}/old_way.trimmed.R2.fq \
                -out_good ${OUTPUT_FOLDER}/old_way.complexity_filtered \
                -out_bad null \
                -lc_method dust \
                -lc_threshold 7 \
                -line_width 0 \
                2>&1

### Renaming Prinseq outputs in order to have somewhat consistent names

for i in 1 2
    do
        mv ${OUTPUT_FOLDER}/old_way.complexity_filtered_${i}.fastq ${OUTPUT_FOLDER}/old_way.complexity_filtered.R${i}.fq
    done

##################################################################################################
############################# bbduk : PhiX decontamination #######################################
##################################################################################################

PHIX_REF=/home/taulagnon/miniconda3/envs/env_test/share/bbmap/resources/phix174_ill.ref.fa.gz

### PhiX filtering for merged (M) reads

echo "Filtering PhiX sequences from" ${OUTPUT_FOLDER}"/old_way.complexity_filtered.M.fq"

bbduk.sh in=${OUTPUT_FOLDER}/old_way.complexity_filtered.M.fastq \
         out=${OUTPUT_FOLDER}/old_way.phiX_filtered.M.fq \
         ref=$PHIX_REF \
         k=31 hdist=1 qin=33 \
         stats=${OUTPUT_FOLDER}/old_way.phiX_removal_stats.M.txt \
         overwrite=t

### PhiX filtering for unmerged (R1+R2) reads

echo "Filtering PhiX sequences from" ${OUTPUT_FOLDER}"/old_way.complexity_filtered.R1.fq and" ${OUTPUT_FOLDER}"/old_way.complexity_filtered.R1.fq"

bbduk.sh in=${OUTPUT_FOLDER}/old_way.complexity_filtered.R1.fq \
         in2=${OUTPUT_FOLDER}/old_way.complexity_filtered.R2.fq \
         out=${OUTPUT_FOLDER}/old_way.phiX_filtered.R1.fq \
         out2=${OUTPUT_FOLDER}/old_way.phiX_filtered.R2.fq \
         ref=$PHIX_REF \
         k=31 hdist=1 qin=33 \
         stats=${OUTPUT_FOLDER}/old_way.phiX_removal_stats.PE.txt \
         overwrite=t

##################################################################################################
############################# Prinseq : Deduplication of reads ###################################
##################################################################################################

### Duplicate Removal for merged (M) reads

echo "Removing duplicates from" ${OUTPUT_FOLDER}"/old_way.phiX_filtered.M.fq"

prinseq-lite.pl -fastq ${OUTPUT_FOLDER}/old_way.phiX_filtered.M.fq \
                -out_format 1 \
                -out_good ${OUTPUT_FOLDER}/old_way.duplicates_removed.M \
                -out_bad null -derep 1 -line_width 0 \
                2> /dev/null

### Duplicate Removal for unmerged (R1+R2) reads

echo "Removing duplicates from" ${OUTPUT_FOLDER}"/old_way.phiX_filtered.R1.fq and" ${OUTPUT_FOLDER}"/old_way.phiX_filtered.R2.fq"

prinseq-lite.pl -fastq ${OUTPUT_FOLDER}/old_way.phiX_filtered.R1.fq \
                -fastq2 ${OUTPUT_FOLDER}/old_way.phiX_filtered.R2.fq \
                -out_format 1 \
                -out_good ${OUTPUT_FOLDER}/old_way.duplicates_removed \
                -out_bad null -derep 1 -line_width 0 \
                2> /dev/null

### Renaming Prinseq outputs in order to have somewhat consistent names

for i in 1 2
    do
        mv ${OUTPUT_FOLDER}/old_way.duplicates_removed_${i}.fasta ${OUTPUT_FOLDER}/old_way.duplicates_removed.R${i}.fa
    done
