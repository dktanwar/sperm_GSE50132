#!/bin/bash	
##This script uses to annotate small RNA
date
echo ""
thread="16"
adapter5="default"
adapter3="CTGTAGGCACCATCAAT"
min_length="15"
max_length="45"
mismatch="1"
input_query_address="/mnt/IM/projects/software/shortRNA_reports/sperm_gapp2014/data/"
input_query_name="SRR955628"
input_query_suffix="fastq"
output_address="/mnt/IM/projects/software/shortRNA_reports/sperm_gapp2014/sports_adapter/1_${input_query_name}/"
script_address="/mnt/IM/DKT/software/downloads/sports1.1/source/"


if [ ! -d "${output_address}" ]; then
	mkdir -p ${output_address}
fi
ln -s ${input_query_address}${input_query_name}.${input_query_suffix} ${output_address}
cd ${output_address}
echo ""
echo "remove 3' adapter"
cutadapt -j ${thread} -a CTGTAGGCACCATCAAT -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
rm -rf ${output_address}${input_query_name}.${input_query_suffix}
		
perl ${script_address}fastq2fasta.pl ${output_address}${input_query_name}_trim_1.${input_query_suffix} > ${output_address}${input_query_name}_trim_1.fa
rm ${output_address}${input_query_name}_trim_1.${input_query_suffix}
input_query_suffix=fa

perl ${script_address}fastaparse.pl ${output_address}${input_query_name}_trim_1.fa -b > ${output_address}${input_query_name}_trim_2.fa 2>${output_address}${input_query_name}_discarded_reads.fa
rm ${output_address}${input_query_name}_trim_1.fa

perl ${script_address}fastaparse.pl ${output_address}${input_query_name}_trim_2.fa -a ${min_length} > ${output_address}${input_query_name}_trim_3.fa 2>${output_address}${input_query_name}_too_short_reads.fa
rm ${output_address}${input_query_name}_trim_2.fa

perl ${script_address}fastaparse.pl ${output_address}${input_query_name}_trim_3.fa -c ${max_length} > ${output_address}${input_query_name}_trim_4.fa 2>${output_address}${input_query_name}_too_long_reads.fa
rm ${output_address}${input_query_name}_trim_3.fa

perl ${script_address}combine_reads.pl ${output_address}${input_query_name}_trim_4.fa > ${output_address}${input_query_name}.fa
rm ${output_address}${input_query_name}_trim_4.fa
		

input=${output_address}${input_query_name}.fa	

###step1: match to genome
echo "match to genome"
bowtie_address=/mnt/IM/reference/genome/gencode/bowtie/mm10
output_match=${output_address}${input_query_name}_match_genome.fa
output_unmatch=${output_address}${input_query_name}_unmatch_genome.fa
output_detail=${output_address}${input_query_name}_output_match_genome

touch ${output_match}
touch ${output_unmatch}
bowtie ${bowtie_address} -f ${input} -v ${mismatch} -k 1 -p ${thread} --al ${output_match} --un ${output_unmatch} > ${output_detail}

input_match=${output_address}${input_query_name}_match_genome.fa
input_unmatch=${output_address}${input_query_name}_unmatch_genome.fa
###step2: match to rRNA database
echo ""
echo "match to rRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_rRNA_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_rRNA_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}
name=rRNA_4.5S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_4.5S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_5S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_5S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_5.8S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_5.8S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_12S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_12S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_16S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_16S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_18S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_18S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_28S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_28S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_45S
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_45S

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_RNY1
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_RNY1

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
name=rRNA_RNY3
bowtie_address=../annotation/Mus_musculus/rRNAdb/mouse_rRNA_RNY3

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
###step3: match to tRNA database
echo ""
echo "match to tRNA database"

######match genome part - tRNA-mature
name=tRNA_mature
bowtie_address=../annotation/Mus_musculus/GtRNAdb/mm10/mm10-tRNAs_CCA
echo ""
echo "match to tRNA_mature-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

input_match=${output_unmatch_match_genome}

######match genome part - tRNA
name=tRNA_pre
bowtie_address=../annotation/Mus_musculus/GtRNAdb/mm10/mm10-tRNAs
echo ""
echo "match to tRNA-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

######unmatch genome part - tRNA-mature
name=tRNA_mature
bowtie_address=../annotation/Mus_musculus/GtRNAdb/mm10/mm10-tRNAs_CCA
echo ""
echo "match to tRNA_mature-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

input_unmatch=${output_unmatch_unmatch_genome}

######unmatch genome part - tRNA
name=tRNA_pre
bowtie_address=../annotation/Mus_musculus/GtRNAdb/mm10/mm10-tRNAs
echo ""
echo "match to tRNA-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
###step4: match to mito_tRNA database
echo ""
echo "match to mito_tRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_mt_tRNA_mature_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_mt_tRNA_mature_unmatch_genome
	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}
name=mt_tRNA_mature
bowtie_address=../annotation/Mus_musculus/GtRNAdb/mm10/mm10-mt_tRNAs_CCA

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
echo ""
output_detail_match_genome=${output_address}${input_query_name}_output_mt_tRNA_pre_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_mt_tRNA_pre_unmatch_genome

if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}
name=mt_tRNA_pre
bowtie_address=../annotation/Mus_musculus/GtRNAdb/mm10/mm10-mt_tRNAs

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
###step5: match to microRNA database
echo ""
echo "match to microRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_miRNA_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_miRNA_unmatch_genome
	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}
name=miRNA
bowtie_address=../annotation/Mus_musculus/miRBase/21/miRBase_21-mmu

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
###step6: match to ensembl database
echo ""
echo "match to ensembl database"
output_detail_match_genome=${output_address}${input_query_name}_output_ensembl_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_ensembl_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}
name=ensembl
bowtie_address=../annotation/Mus_musculus/Ensembl/release-89/Mus_musculus.GRCm38.ncrna

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
###step7: match to rfam database
echo ""
echo "match to rfam database"
output_detail_match_genome=${output_address}${input_query_name}_output_rfam_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_rfam_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}
name=rfam
bowtie_address=../annotation/Mus_musculus/Rfam/12.3/Rfam-12.3-mouse

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
###step8: match to piRNA database
echo ""
echo "match to piRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_piRNA_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_piRNA_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}
name=piRNA
bowtie_address=../annotation/Mus_musculus/piRBase/piR_mouse

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}
perl ${script_address}annotation.pl ${output_address}${input_query_name}

if [ ! -d "${output_address}${input_query_name}_result" ]; then
	mkdir -p ${output_address}${input_query_name}_result
fi

if [ ! -d "${output_address}${input_query_name}_processed" ]; then
	mkdir -p ${output_address}${input_query_name}_processed
fi

if [ ! -d "${output_address}${input_query_name}_fa" ]; then
	mkdir -p ${output_address}${input_query_name}_fa
fi

mv ${output_address}${input_query_name}_*.txt ${output_address}${input_query_name}_result
mv ${output_address}${input_query_name}_output_* ${output_address}${input_query_name}_processed
mv ${output_address}${input_query_name}*.f* ${output_address}${input_query_name}_fa
echo ""
echo "generating graph"
echo ""
######overall length distribution figure
Rscript --vanilla ${script_address}overall_RNA_length_distribution.R ${output_address} ${input_query_name}
echo "overall length distribution figure generated"

######pre-tRNA mapping figure
cat ${output_address}${input_query_name}_processed/${input_query_name}_output_*tRNA_pre_*_genome > ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_pre
if [ -s ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_pre ]; then
	perl ${script_address}tRNA_mapping.pl ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_pre ${output_address}${input_query_name}_result/${input_query_name}_summary.txt > ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_pre_mapping.txt
	Rscript --vanilla ${script_address}tRNA_mapping.R ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_pre_mapping.txt ${output_address}${input_query_name}_result/${input_query_name}_tRNA_pre_mapping.pdf
	echo "pre-tRNA mapping figure generated"
fi

######mature-tRNA mapping figure
cat ${output_address}${input_query_name}_processed/${input_query_name}_output_*tRNA_mature_*_genome > ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_mature
if [ -s ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_mature ]; then
	perl ${script_address}tRNA_mapping.pl ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_mature ${output_address}${input_query_name}_result/${input_query_name}_summary.txt > ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_mature_mapping.txt
	Rscript --vanilla ${script_address}tRNA_mapping.R ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_mature_mapping.txt ${output_address}${input_query_name}_result/${input_query_name}_tRNA_mature_mapping.pdf
	echo "mature-tRNA mapping figure generated"
fi
temp_length=4.5S=174,5S=121,5.8S=157,12S=954,16S=1583,18S=1870,28S=4730,45S=13400,RNY1=112,RNY3=102
######rRNA length distribution figure
Rscript --vanilla ${script_address}rRNA_length_distribution.R ${output_address} ${input_query_name} ${temp_length}
echo "rRNA length distribution figure generated"

######rRNA mapping figure
Rscript --vanilla ${script_address}rRNA_mapping.R ${output_address} ${input_query_name} ${temp_length}
echo "rRNA mapping figure generated"
rm ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA
rm ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_mapping.txt
rm -rf ${output_address}${input_query_name}_fa
rm -rf ${output_address}${input_query_name}_processed
echo ""
date