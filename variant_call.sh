#!/bin/bash
set -eo pipefail

# Variable definitions
reference=/path/to/reference/file/$1
prefix=/path/to/folder/$2 
read1fn="/path/to/trimmed/sequences/$2_R1_trimmed.fastq"
read2fn="/path/to/trimmed/sequences/$2_R2_trimmed.fastq"

cd $prefix
if [ ! -f ../*.dict ]
then gatk CreateSequenceDictionary -R $reference
fi 

bwa index $reference
samtools faidx $reference

# Aligning the trimmed sequences to the reference and mapping
bwa mem $reference $read1fn $read2fn | samtools view -bS - | samtools sort - -o "$prefix.sorted.bam"
gatk FastqToSam -F1 $read1fn -F2 $read2fn -O $prefix.unmapped.bam -SM $prefix.sorted.bam

# Replacing read groups with mapped and unmapped BAM files using prepared libraries and sequencing information
gatk AddOrReplaceReadGroups -I  $prefix.sorted.bam -O $prefix.sorted-RG.bam -RGID 2 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $prefix
gatk AddOrReplaceReadGroups -I  $prefix.unmapped.bam -O $prefix.unmapped-RG.bam -RGID 2 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $prefix

# Merging mapped and unmapped BAM files
gatk MergeBamAlignment --ALIGNED_BAM $prefix.sorted-RG.bam --UNMAPPED_BAM $prefix.unmapped-RG.bam -O $prefix.merged.bam -R $reference

# Removing duplicate sequences
gatk MarkDuplicates -I $prefix.merged.bam -O $prefix.marked.bam -M $prefix.metrics.txt
samtools index $prefix.marked.bam

# Creating the GVCF file 
gatk HaplotypeCaller -I $prefix.marked.bam -O $prefix-g.vcf -ERC GVCF -R $reference

# Removing intermediate files
rm $prefix.sorted.bam $prefix.unmapped.bam $prefix.merged.bam $prefix.unmapped-RG.bam $prefix.sorted-RG.bam
