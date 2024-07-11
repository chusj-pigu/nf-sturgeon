params.simplex = true
params.dorado_cpu = false
params.b = null

include { basecall } from './subworkflows/dorado'
include { qs_filter } from './subworkflows/samtools'
include { ubam_to_fastq as ubam_to_fastq_p } from './subworkflows/samtools'
include { ubam_to_fastq as ubam_to_fastq_f } from './subworkflows/samtools'
include { nanoplot } from "./subworkflows/nanoplot"
include { mapping as map_hg38 } from './subworkflows/minimap'
include { mapping as map_chm13 } from './subworkflows/minimap'
include { sam_sort } from './subworkflows/samtools'
include { mosdepth } from './subworkflows/mosdepth'
include { multiqc } from './subworkflows/multiqc'
include { gather_sturgeon } from './subworkflows/ingress'
include { adjust_mods } from './subworkflows/modkit'
include { extract } from './subworkflows/modkit'
include { predict } from './subworkflows/sturgeon'
include { inputtobed } from './subworkflows/sturgeon'

workflow {
    pod5_ch = Channel.fromPath(params.pod5)
    ref_hg38_ch = Channel.fromPath(params.ref_hg38)
    ref_hgchm13_ch = Channel.fromPath(params.ref_hgchm13)
    d_model = Channel.fromPath(params.dorado_model)
    stmodel_ch = Channel.fromPath(params.sturgeon_model)

    basecall(pod5_ch, d_model)

    qs_filter(basecall.out)
    nanoplot(basecall.out)

    fq_pass = ubam_to_fastq_p(qs_filter.out.ubam_pass)
    fq_fail = ubam_to_fastq_f(qs_filter.out.ubam_fail)

    map_hg38(ref_hg38_ch, fq_pass)
    map_chm13(ref_hgchm13_ch, fq_pass)

    map_ch = Channel.fromList(map_hg38.out, map_chm13.out)

    sam_sort(map_ch)
    mosdepth(sam_sort.out)

    multi_ch = Channel.empty()
        .mix(nanoplot.out,mosdepth.out)
        .collect()
    multiqc(multi_ch)

    ch13_sorted_ch = sam_sort.out
        .map { hg38, ch13 -> ch13 }
        .view()

   // adjust_mod(ch13_sorted_ch)
    // extract(adjust_mod.out)
   // gather(adjust_mod.out, extract.out)
   // sturgeon_bed(gather.out)
   // sturgeon_predict(sturgeon_bed.out, model_ch)
}
