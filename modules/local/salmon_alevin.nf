
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SALMON_ALEVIN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::salmon=1.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/salmon:1.4.0--h84f40af_1"
    } else {
        container "quay.io/biocontainers/salmon:1.4.0--h84f40af_1"
    }

    input:
    tuple val(meta), path(reads)
    path index
    path txp2gene
    val protocol
    path whitelist

    output:
    tuple val(meta), path("*_alevin_results"), emit: alevin_results
    path "*.version.txt"                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    salmon alevin \\
        -l ISR \\
        -p $task.cpus \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --${protocol} \\
        -i $index \\
        --tgMap $txp2gene \\
        --dumpFeatures --dumpMtx \\
        $options.args \\
        -o ${prefix}_alevin_results

    mv ${whitelist} ${prefix}_alevin_results/alevin/whitelist.txt

    salmon --version | sed -e "s/salmon //g" > ${software}.version.txt
    """
}