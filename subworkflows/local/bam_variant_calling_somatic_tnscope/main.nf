//
//
// SENTIEON TNSCOPE: tumor-normal mode variantcalling
//

include { SENTIEON_TNSCOPE } from '../../../modules/nf-core/sentieon/tnscope/main' 
include { GATK4_MERGEVCFS  as MERGE_TNSCOPE } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_SOMATIC_TNSCOPE {
    take:
    input                     // channel: [ meta, [ input ], [ input_index ] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                     // channel: /path/to/reference/dict
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index
    intervals                 // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for TNSCOPE module
        .map{ meta, input, index, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input, index, intervals ] }
        .view()

    //If no germline resource is provided, then create an empty channel to avoid GetPileupsummaries from being run
    //germline_resource_pileup     = (germline_resource && germline_resource_tbi) ? germline_resource : Channel.empty()
    //germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()

    // Separate normal cram files
    /* Extract tumor cram files
    ch_cram = input.multiMap{ meta, cram, crai ->
            normal: [ meta - meta.subMap('tumor_id') , cram[0], crai[0] ]
            tumor:  [ meta - meta.subMap('tumor_id') , cram[1], crai[1] ]
        }

    // Remove duplicates from normal channel and merge normal and tumor crams by patient
    ch_tn_cram =  ch_cram.normal.unique().mix(ch_cram.tumor).groupTuple()
    */
    SENTIEON_TNSCOPE( 
        input_intervals,
        fasta, 
        fai, 
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        [], //cosmic
        [] // cosmic_tbi
    )
    versions = versions.mix(SENTIEON_TNSCOPE.out.versions)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = SENTIEON_TNSCOPE.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }
    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = SENTIEON_TNSCOPE.out.index.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple()

    // Merge if required
    MERGE_TNSCOPE(vcf_to_merge, dict)
    versions = versions.mix(MERGE_TNSCOPE.out.versions)
    
    // Mix intervals and no_intervals channels together
    // Remove unnecessary metadata and add variantcaller
    vcf   = Channel.empty()
        .mix(MERGE_TNSCOPE.out.vcf, vcf_branch.no_intervals)
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'sentieon_tnscope' ], vcf ] }

    index = Channel.empty()
        .mix(MERGE_TNSCOPE.out.tbi, tbi_branch.no_intervals)
        .map{ meta, tbi -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'sentieon_tnscope' ], tbi ] }

    emit:
    vcf      // channel: [ meta, vcf ]
    index    // channel: [ meta, index ]
    versions // channel: [ versions.yml ]
}
