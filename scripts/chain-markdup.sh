#!/usr/bin/env bash
set -vexubE -o pipefail

# check input must be (aligned.bam, template.bam, out.bam)
if [ $# != 3 ] 
then
    echo 'chain-markdup.sh aligned.bam template_or_subreads.bam output.bam'
    echo 'e.g., chain-markdup.sh ngmlr.bam subreads.bam output.bam' 
    exit 0
fi

input_bam=$1
template_bam=$2
output_bam=$3

# input_bam=mm2.bam
# output_bam=out.bam
# template_bam=tiny.subreads.bam

tmp_bam=${output_bam}.markdup.in.bam
markdup_bam=${output_bam}.markdup.out.bam
chain_sam=${output_bam}.chain.sam

# mkdir -p msrt && (samtools view -H ${input_bam} ; (samtools view ${input_bam} | sort --stable -T msrt -k1,1)) | samtools sort --threads 16 - > ${tmp_bam} 
pbsvutil markduplicates ${input_bam} ${markdup_bam} --template_bam ${template_bam}
mkdir -p csrt ; samtools view ${markdup_bam} | sort --stable -k1,1 | cat <(samtools view -H ${markdup_bam} -) - | pbsvutil chain - ${chain_sam}
samtools sort --threads 4  ${chain_sam} -o ${output_bam} ; samtools index ${output_bam}
