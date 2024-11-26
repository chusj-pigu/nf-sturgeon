include { basecall } from '../modules/dorado'
include { qs_filter } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_p } from '../modules/samtools'
include { ubam_to_fastq as ubam_to_fastq_f } from '../modules/samtools'
include { nanoplot } from "../modules/nanoplot"
include { ALIGNMENT } from '../subworkflows/mapping'
include { multiqc } from '../modules/multiqc'

workflow SIMPLEX {
    take:
    pod5
    ubam
    model
    
    main:
    basecall(pod5, model, ubam)

    qs_filter(basecall.out)
    nanoplot(basecall.out)

    ubam_to_fastq_p(qs_filter.out.ubam_pass)
    ubam_to_fastq_f(qs_filter.out.ubam_fail)

    emit:
    fq_pass = ubam_to_fastq_p.out
    nanoplot_res = nanoplot.out

}