#!/usr/bin/bash

CELLRANGER_PATH="/path_to/cellranger/cellranger-2.1.1"
transcriptome="/path_to/cellranger/refdata/refdata-cellranger-mm10-1.2.0"

preamble="PATH=$CELLRANGER_PATH:\$PATH ; export PATH ; cd \$PBS_O_WORKDIR"

TIME="48:00:00"
THREADS="16"
MEM="150gb"

if [[ $# -lt "5" ]]
then
	echo "Usage: bash run_CellRanger_count.sh <fastq_files> <id> <project> <sample_names> <cells> [<f>]"
	echo "<fastq_files>      output path prefix from \"cellranger mkfastq\" run, one per line"
	echo "<id>               unique ID for this CellRanger run"
	echo "<project>          name of the project"
	echo "<sample_name>      should be found inside <fastq_files> (the same used in mkfastq!!!)"
	echo "<cells>            number of expected cells"
	echo "<f>                force pipeline to use this number of cells, bypassing the cell detection algorithm [OPTIONAL - set \"f\" to enable]"
	exit
fi

fasta_files=$1
id=$2
project=$3
sample_name=$4
cells=$5
force="n"
if [[ $# -gt "5" ]]
then
	force=$6
fi

fastqs=$(cat $fasta_files | tr "\n" "," | sed "s/,$//g")

if [[ $force != "f" ]]
then
	command="cellranger count --id $id --project $project --fastqs $fastqs --transcriptome $transcriptome --sample $sample_name --expect-cells $cells > count_$id.STDOUT 2> count_$id.STDERR"
else 
	command="cellranger count --id $id --project $project --fastqs $fastqs --transcriptome $transcriptome --sample $sample_name --force-cells $cells > count_$id.STDOUT 2> count_$id.STDERR"
fi
echo "$preamble ; $command"
echo "$preamble ; $command" | qsub -q batch -N cellranger_count_$id -l walltime=$TIME,nodes=1:ppn=$THREADS,mem=$MEM

exit
