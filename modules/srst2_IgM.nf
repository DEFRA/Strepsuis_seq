process SRST2_IGM {
    tag "$pairId"
    errorStrategy 'ignore'
    maxForks 2
    publishDir "${params.outdir}/Results_${params.dataDir}_${params.today}/protease", mode: "copy"

    input:
        tuple val(pairId), file(R1), file(R2)

    output:
	    tuple val(pairId), file("IgMprotease_${pairId}_*.txt"), emit: protease_results
    
    script:
    """
    srst2.py --input_pe $R1 $R2 --forward _.*_R1 --reverse _.*_R2 \
             --output IgMprotease_$pairId --mlst_db ${params.igm_db} --mlst_definitions ${params.igm_def}\
             --mlst_delimiter "-" --max_unaligned_overlap 75 --min_coverage 60 --max_divergence 20 --log
    """
}