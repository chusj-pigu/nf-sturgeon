params.dorado_cpu = false
params.b = null
params.no_mod = false

// Usage help

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/nf-sturgeon --pod5 path/to/pod5 --sample_ID "ID"

        Mandatory arguments:
         --pod5                         Path to the directory containing pod5 files
         --ref_hg38                     Path to the GRCh38 genome assembly. If using -profile drac, you can omit this option and the reference under Platform directory in ctb-noncodo will be used by default.
         --ref_hgchm13                  Path to the Telomere-to-telomere reference genome (CHM13v2). If using -profile drac, you can omit this option and the reference under Platform directory in ctb-noncodo will be used by default.
         --sturgeon_model               Path to the sturgeon model. If using -profile drac, you can omit this option and the reference under Platform directory in ctb-noncodo will be used by default.

         Optional arguments:
         --dorado_model                 Basecalling model, path is required when running with drac profile [default: path to sup@v5.0.0]
         --out_dir                      Output directory to place mapped files and reports in [default: output]
         --sample_id                    Will name output files according to sample id [default: reads]
         --m_bases                      Modified bases to be called, separated by commas if more than one is desired. Requires path to model if run with drac profile [default: 5mCG_5hmCG].
         --m_bases_path                 Path for the modified basecalling model, required when running with drac profile [default: path to sup@v5.0.0_5mCG_5hmCG]
         -profile                       Use standard for running locally, or drac when running on Digital Research Alliance of Canada Narval [default: standard]
         --threads                      Number of threads to use for mapping [default: 40]
         --skip_basecall                Will skip basecalling, takes fastq as input. 
         --skip_hg38                    Will skip parallel mapping to hg38.
         --help                         This usage statement.
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

include { basecall } from './subworkflows/dorado'
include { qs_filter } from './subworkflows/samtools'
include { ubam_to_fastq as ubam_to_fastq_p } from './subworkflows/samtools'
include { ubam_to_fastq as ubam_to_fastq_f } from './subworkflows/samtools'
include { nanoplot } from "./subworkflows/nanoplot"
include { mapping as map_hg38 } from './subworkflows/minimap'
include { mapping as map_chm13 } from './subworkflows/minimap'
include { sam_sort as sort_hg38 } from './subworkflows/samtools'
include { sam_sort as sort_chm13 } from './subworkflows/samtools'
include { mosdepth as mos_hg38 } from './subworkflows/mosdepth'
include { mosdepth as mos_chm13 } from './subworkflows/mosdepth'
include { multiqc } from './subworkflows/multiqc'
include { gather_sturgeon } from './subworkflows/ingress'
include { adjust_mods } from './subworkflows/modkit'
include { extract } from './subworkflows/modkit'
include { predict } from './subworkflows/sturgeon'
include { inputtobed } from './subworkflows/sturgeon'

workflow {

    ref_hgchm13_ch = Channel.fromPath(params.ref_hgchm13)
    d_model = Channel.fromPath(params.dorado_model)
    stmodel_ch = Channel.fromPath(params.sturgeon_model)

    if (params.skip_basecall && params.skip_hg38) {
        fq_pass = Channel.fromPath(params.fastq)
        
        map_chm13(ref_hgchm13_ch, fq_pass)
        sort_chm13(map_chm13.out)
        mos_chm13(sort_chm13.out)

        multi_ch = Channel.empty()
            .mix(mos_chm13.out)
            .collect()
        multiqc(multi_ch)
    } else if (!params.skip_basecall && params.skip_hg38) {
        pod5_ch = Channel.fromPath(params.pod5)
        basecall(pod5_ch, d_model)

        qs_filter(basecall.out)
        nanoplot(basecall.out)

        fq_pass = ubam_to_fastq_p(qs_filter.out.ubam_pass)
        fq_fail = ubam_to_fastq_f(qs_filter.out.ubam_fail)

        map_chm13(ref_hgchm13_ch, fq_pass)
        sort_chm13(map_chm13.out)
        mos_chm13(sort_chm13.out)

        multi_ch = Channel.empty()
            .mix(nanoplot.out,mos_chm13.out)
            .collect()
        multiqc(multi_ch)
    } else if (!params.skip_hg38 && params.skip_basecall) {
        fq_pass = Channel.fromPath(params.fastq)

        map_hg38(ref_hg38_ch, fq_pass)
        map_chm13(ref_hgchm13_ch, fq_pass)

        sort_hg38(map_hg38.out)
        sort_chm13(map_chm13.out)

        mos_hg38(sort_hg38.out)
        mos_chm13(sort_chm13.out)

        multi_ch = Channel.empty()
            .mix(mos_hg38.out,mos_chm13.out)
            .collect()
        multiqc(multi_ch)

    } else {
        pod5_ch = Channel.fromPath(params.pod5)
        basecall(pod5_ch, d_model)

        qs_filter(basecall.out)
        nanoplot(basecall.out)

        fq_pass = ubam_to_fastq_p(qs_filter.out.ubam_pass)
        fq_fail = ubam_to_fastq_f(qs_filter.out.ubam_fail)

        map_hg38(ref_hg38_ch, fq_pass)
        map_chm13(ref_hgchm13_ch, fq_pass)

        sort_hg38(map_hg38.out)
        sort_chm13(map_chm13.out)

        mos_hg38(sort_hg38.out)
        mos_chm13(sort_chm13.out)

        multi_ch = Channel.empty()
            .mix(nanoplot.out,mos_hg38.out,mos_chm13.out)
            .collect()
        multiqc(multi_ch)

    }

    bam_only_chm13 = sort_chm13.out
        .map { bam, bai -> bam }

    adjust_mods(bam_only_chm13)
    extract(adjust_mods.out)
    gather_sturgeon(adjust_mods.out, extract.out)
    inputtobed(gather_sturgeon.out)
    predict(inputtobed.out, stmodel_ch)
}
