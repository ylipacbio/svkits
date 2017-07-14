#!/usr/bin/env bash
: ${1?"Usage: start end"}
: ${2?"Usage: start end"}

ref_fa=/home/UNIXHOME/yli/projects/svkits/ecoli/sim_reference/action_delete_50_id_0/reference.fasta
sr_bam=/pbi/collections/315/3150572/r54006_20170110_165905/2_B01/m54006_170111_014001.subreads.bam
in_bam=/home/UNIXHOME/yli/projects/svkits/ecoli/sim_reference/action_delete_50_id_0/pbsmrtpipe/tasks/pbsvtools.tasks.gather_align-1/alignments.bam

start=$1
end=$2
echo start=$start, end=$end

o_prefix=reads_span_${start}_${end}
subset_org_sam=${o_prefix}.original.sam
subset_org_bam=${o_prefix}.original.bam
subset_zmws=${o_prefix}.zmws

samtools index $in_bam

# samtools view region
samtools view -H $in_bam > ${subset_org_sam}
samtools view $in_bam "000000F|arrow:${start}-${end}" >> ${subset_org_sam}
samtools view -bS ${subset_org_sam} -o ${subset_org_bam}
samtools index ${subset_org_bam}
cat ${subset_org_sam} |grep -v '@' | cut -f 1 > ${subset_zmws}

subset_sr_bam=$o_prefix.subreads.bam
blasr_out=${o_prefix}.blasr.bam
blasr_out_m5=${o_prefix}.blasr.m5
ngmlr_out=${o_prefix}.ngmlr.bam

make-subreads-bam-from-zmws $subset_zmws $sr_bam $subset_sr_bam
blasr $subset_sr_bam $ref_fa --out $blasr_out --bam 
blasr $subset_sr_bam $ref_fa --out $blasr_out_m5 -m 5
samtools sort $blasr_out > ${blasr_out}.tmp && mv ${blasr_out}.tmp  $blasr_out
samtools index $blasr_out
pbsv align $ref_fa $subset_sr_bam $ngmlr_out
