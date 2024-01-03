process SIMPLEAF_QUANT_AUTO {
    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::simpleaf=0.10.0-1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.10.0--h9f5acd7_1' :
        'biocontainers/simpleaf:0.10.0--h9f5acd7_1' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path index
    path txp2gene
    path whitelist

    output:
    tuple val(meta), path("*_alevin_results"), emit: alevin_results
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def args_list = args.tokenize()
    def prefix    = task.ext.prefix ?: "${meta.id}"

    unfiltered_command = ""
    save_whitelist     = ""
    if (!(args_list.any { it in ['-k', '--knee', '-e', '--expect-cells', '-f', '--forced-cells']} || meta.expected_cells)) {
        unfiltered_command = "-u whitelist.uncompressed.txt"
        save_whitelist     = "mv whitelist.uncompressed.txt ${prefix}_alevin_results/"
    }

    // expected cells
    def expect_cells = meta.expected_cells ? "--expect-cells $meta.expected_cells" : ''

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()

    """

    # extract barcode whitelists for three different versions (V1, V2, V3).
    for PROTOCOL in 1 2 3; do
        gzip -dcf "whitelist/10x_V\${PROTOCOL}_barcode_whitelist.txt.gz" |
        sed "s/\$/ \$PROTOCOL/" # append column containing protocol version
    done > barcodes

    # detect chemistry version (1/2/3)
    PROTOCOL=\$(
        gzip -dcf "${forward[0]}" |
        awk 'FNR % 4 == 2' | # extract the read sequence from FastQ
        head -n 100000 | # the first 100k reads should suffice
        awk '
            { \$1 = substr(\$1, 1, 14) } # the barcode is in the first 14 bases (we ignore that V2/3 chemistries have longer barcodes)
            FILENAME ~ /barcodes/ { protocols[\$2]; barcodes[\$0] } # cache barcodes in memory
            FILENAME ~ /stdin/ { for (p in protocols) if (\$0" "p in barcodes) count[p]++; total++ } # count matches for each chemistry
            END { for (p in protocols) if (count[p]/total > 0.85) print p } # output chemistries with >85% matches
        ' barcodes /dev/stdin || true
    )

    # based on the chemistry to choose the correct 10x barcode whitelist
    gzip -dcf "whitelist/10x_V\${PROTOCOL}_barcode_whitelist.txt.gz" > whitelist.uncompressed.txt

    # export required var
    export ALEVIN_FRY_HOME=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf quant
    simpleaf quant \\
        -1 ${forward.join( "," )} \\
        -2 ${reverse.join( "," )} \\
        -i ${index} \\
        -o ${prefix}_alevin_results \\
        -m $txp2gene \\
        -t $task.cpus \\
        -c "10xv\$PROTOCOL" \\
        $expect_cells \\
        $unfiltered_command \\
        $args

    $save_whitelist
    [[ ! -f ${prefix}_alevin_results/af_quant/all_freq.bin ]] && cp ${prefix}_alevin_results/af_quant/permit_freq.bin ${prefix}_alevin_results/af_quant/all_freq.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
