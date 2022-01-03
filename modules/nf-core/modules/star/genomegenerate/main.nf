// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::star=2.7.8a' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/star:2.7.8a--h9ee0642_1'
    } else {
        container 'quay.io/biocontainers/star:2.7.8a--h9ee0642_1'
    }
    
    input:
    path fasta
    path gtf

    output:
    path "star"         , emit: index
    path "*.version.txt", emit: version

    script:
    def software  = getSoftwareName(task.process)
    def memory    = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --runThreadN $task.cpus \\
        $memory \\
        $options.args

    STAR --version | sed -e "s/STAR_//g" > ${software}.version.txt
    """
}
