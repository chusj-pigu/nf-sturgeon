include { gather_sturgeon } from '../modules/ingress'
include { adjust_mods } from '../modules/modkit'
include { extract } from '../modules/modkit'
include { predict } from '../modules/sturgeon'
include { inputtobed } from '../modules/sturgeon'

workflow CNS_RESULT {
    take:
    bam_chm13
    sturgeon_model

    main:
    adjust_mods(bam_chm13)
    extract(adjust_mods.out)
    gather_sturgeon(adjust_mods.out, extract.out)
    inputtobed(gather_sturgeon.out)
    predict(inputtobed.out, sturgeon_model)
}