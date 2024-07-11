module load nextflow
module load apptainer

BAM="/project/ctb-noncodo/Platform/2024/05/PPT45/sturgeon/SIGN1048.chm13_aligned_sorted.bam"
NAME=SIGN1048

nextflow run main.nf -resume \
    --bam $BAM \
    --sample_id $NAME