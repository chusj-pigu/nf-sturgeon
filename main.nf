
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

include { SIMPLEX } from "./subworkflows/simplex"
include { ALIGNMENT as ALIGN_hg38 } from './subworkflows/mapping'
include { ALIGNMENT as ALIGN_chm13 } from './subworkflows/mapping'
include { multiqc } from './modules/multiqc'
include { CNS_RESULT } from './subworkflows/sturgeon'

workflow {

    ref_hgchm13_ch = Channel.fromPath(params.ref_hgchm13)
    d_model = Channel.fromPath(params.dorado_model)
    stmodel_ch = Channel.fromPath(params.sturgeon_model)

    if (params.skip_basecall && params.skip_hg38) {
        fq_pass = Channel.fromPath(params.fastq)
        
        ALIGN_chm13(fq_pass, ref_hgchm13_ch)

        multi_ch = Channel.empty()
            .mix(ALIGN_chm13.out.mosdepth_dist, ALIGN_chm13.out.mosdepth_summary, ALIGN_chm13.out.mosdepth_bed)
            .collect()
        multiqc(multi_ch)
        
    } else if (!params.skip_basecall && params.skip_hg38) {
        pod5_ch = Channel.fromPath(params.pod5)
        SIMPLEX(pod5_ch,d_model)
        ALIGN_chm13(SIMPLEX.out.fq_pass,ref_hgchm13_ch)

        multi_ch = Channel.empty()
            .mix(SIMPLEX.out.nanoplot_res, ALIGN_chm13.out.mosdepth_dist, ALIGN_chm13.out.mosdepth_summary, ALIGN_chm13.out.mosdepth_bed)
            .collect()
        multiqc(multi_ch)

    } else if (!params.skip_hg38 && params.skip_basecall) {
        fq_pass = Channel.fromPath(params.fastq)
        ref_hg38_ch = Channel.fromPath(params.ref_hg38)

        ALIGN_chm13(fq_pass, ref_hgchm13_ch)
        ALIGN_hg38(fq_pass, ref_hg38_ch)

        multi_ch = Channel.empty()
            .mix(ALIGN_chm13.out.mosdepth_dist, ALIGN_chm13.out.mosdepth_summary, ALIGN_chm13.out.mosdepth_bed, ALIGN_hg38.out.mosdepth_dist, ALIGN_hg38.out.mosdepth_summary, ALIGN_hg38.out.mosdepth_bed)
            .collect()
        multiqc(multi_ch)

    } else {
        pod5_ch = Channel.fromPath(params.pod5)
        ref_hg38_ch = Channel.fromPath(params.ref_hg38)

        SIMPLEX(pod5_ch,d_model)

        ALIGN_chm13(SIMPLEX.out.fq_pass, ref_hgchm13_ch)
        ALIGN_hg38(SIMPLEX.out.fq_pass, ref_hg38_ch)

        multi_ch = Channel.empty()
            .mix(SIMPLEX.out.nanoplot_res,ALIGN_chm13.out.mosdepth_dist, ALIGN_chm13.out.mosdepth_summary, ALIGN_chm13.out.mosdepth_bed, ALIGN_hg38.out.mosdepth_dist, ALIGN_hg38.out.mosdepth_summary, ALIGN_hg38.out.mosdepth_bed)
            .collect()
        multiqc(multi_ch)

    }

    bam_only_chm13 = ALIGN_chm13.out.bam
        .map { bam, bai -> bam }

    CNS_RESULT(bam_only_chm13, stmodel_ch)
}
