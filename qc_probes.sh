#!/bin/bash

module load blat



###CHANGE THE INPUT FILES###
PROBES="../out-filtered-baits.fa"
TRANS="../matricaria_trans.fasta"
CONTIG="../a_radiatus_contigs.fasta"
PLASTID="../a_radiatus_plastid.fasta"
MITO="../helianthus_mitochondrion.fasta"
RDNA="../a_valentinus_rdna.fasta"


INITIAL_PROBES=$(cat $PROBES | grep -c ">")


echo -e "Checking for multiple hit of the probe set\n" 

echo -e "\n\tAligning probes against the transcriptome file..."
blat -t=dna -q=dna -minIdentity=90 -noHead -tileSize=15 -out=pslx  $PROBES $TRANS probes_v_trans.pslx
echo -e "\tFiltering the results..."
cat  probes_v_trans.pslx | awk '{if ($1>110) print $14}' |sort | uniq -c |sort | awk '{if ($1>1)  print $2}'> probes_with_multiple_hit_trans.list
PROBES_DUPLICATED=$(wc -l probes_with_multiple_hit_trans.list|awk '{print $1}')
echo -e "\t At least $PROBES_DUPLICATED probes have duplicated hits, they should be discarded.\n\tList of probes to be discarded in probes_with_multiple_hit_trans.list"


echo -e "\n\tAligning probes against the contig file..."
blat -t=dna -q=dna -minIdentity=90 -noHead -tileSize=15 -out=pslx  $PROBES $CONTIG probes_v_contigs.pslx
echo -e "\tFiltering the results..."
cat  probes_v_contigs.pslx | awk '{if ($1>110) print $14}' |sort | uniq -c |sort | awk '{if ($1>1)  print $2}'> probes_with_multiple_hit_contig.list
PROBES_DUPLICATED=$(wc -l probes_with_multiple_hit_contig.list|awk '{print $1}')
echo -e "\t At least $PROBES_DUPLICATED probes have duplicated hits, they should be discarded.\n\tList of probes to be discarded in probes_with_multiple_hit_contig.list"

echo -e "\n\tAligning probes against the plastid genome..."
blat -t=dna -q=dna -minIdentity=90 -noHead -tileSize=15 -out=pslx  $PROBES $PLASTID probes_v_plastid.pslx
echo -e "\tFiltering the results..."
cat  probes_v_plastid.pslx | awk '{if ($1>110) print $14}' |sort | uniq -c |sort | awk '{if ($1>1)  print $2}' > probes_with_multiple_hit_plastid.list
PROBES_DUPLICATED=$(wc -l probes_with_multiple_hit_plastid.list|awk '{print $1}')
echo -e "\t  $PROBES_DUPLICATED probes might recover plastid reads, they should be discarded.\n\tList of probes to be discarded in probes_with_multiple_hit_plastid.list"

echo -e "\n\tAligning probes against the mitchondrial genome..."
blat -t=dna -q=dna -minIdentity=90 -noHead -tileSize=15 -out=pslx  $PROBES $MITO probes_v_mito.pslx
echo -e "\tFiltering the results..."
cat  probes_v_mito.pslx | awk '{if ($1>110) print $14}' |sort | uniq -c |sort | awk '{if ($1>1)  print $2}' > probes_with_multiple_hit_mitchondrial.list
PROBES_DUPLICATED=$(wc -l probes_with_multiple_hit_mitchondrial.list|awk '{print $1}')
echo -e "\t $PROBES_DUPLICATED probes might recover mitochondrial reads, they should be discarded.\n\tList of probes to be discarded in probes_with_multiple_hit_mitchondrial.list"

echo -e "\n\tAligning probes against the rDNA ..."
blat -t=dna -q=dna -minIdentity=90 -noHead -tileSize=15 -out=pslx  $PROBES $RDNA probes_v_rdna.pslx
echo -e "\tFiltering the results..."
cat  probes_v_rdna.pslx | awk '{if ($1>110) print $14}' |sort | uniq -c |sort | awk '{if ($1>1)  print $2}' > probes_with_multiple_hit_rDNA.list
PROBES_DUPLICATED=$(wc -l probes_with_multiple_hit_rDNA.list|awk '{print $1}')
echo -e "\t $PROBES_DUPLICATED probes might recover rDNA reads, they should be discarded.\n\tList of probes to be discarded in probes_with_multiple_hit_rDNA.list"

cat probes_with_multiple*.list | sort | uniq > final_list_discarded_probes.list
LENGTH=$(wc -l final_list_discarded_probes.list |awk '{print $1}')
echo -e "\n\nFinal list of probes to be discarded in final_list_of_probes.list"
echo "On the $INITIAL_PROBES initial probes set, $LENGTH probes should be discarded"

echo -e "\n\nFiltering markers... in list_markers_with_duplicated_probes.tab"
cat probes_v_contigs.pslx |  awk '{if ($1>110) print $14}' | sort | uniq -c | sort | awk '{if ($1>1) print $2}' | tr "_" "\t" | cut -f 1 | sort | uniq -c | sort -n > list_markers_with_duplicated_probes.tab

