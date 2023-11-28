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
    val protocol
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
    
    //
    // check if users are using one of the mutually excludable parameters:
    //    e.g -k,--knee | -e,--expect-cells | -f, --forced-cells
    //
    unzip_whitelist = ""
    unfiltered_command = ""
    save_whitelist     = ""
    if (!(args_list.any { it in ['-k', '--knee', '-e', '--expect-cells', '-f', '--forced-cells']} || meta.expected_cells)) {
        if (whitelist) {
            //unzip_whitelist = "gzip -dcf $whitelist > whitelist.uncompressed.txt"
            unfiltered_command = "-u whitelist.uncompressed.txt"
            save_whitelist     = "mv whitelist.uncompressed.txt ${prefix}_alevin_results/"
        }
    }


    // expected cells
    def expect_cells = meta.expected_cells ? "--expect-cells $meta.expected_cells" : ''

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()

    """
    # Check if there are more than one fastq files as input
    IFS=','
    read -ra files <<< ${forward.join( "," )}
    length=\${#files[@]}
    
    if [ \${#files[@]} -ge 2 ]; then
        result="\${files[0]}"
    else
        result=${forward.join( "," )}
    fi
    echo "\$result"

    # Downloads and extracts Barcode Whitelists for three different versions (V1, V2, V3).

    
    V1_BARCODES=\$(gunzip -cd \$PWD/whitelist/10x_V1_barcode_whitelist.txt.gz)
    V2_BARCODES=\$(gunzip -cd \$PWD/whitelist/10x_V2_barcode_whitelist.txt.gz | cut -c1-14)
    V3_BARCODES=\$(gunzip -cd \$PWD/whitelist/10x_V3_barcode_whitelist.txt.gz | cut -c1-14)

    

    # Output chemistry version (V1, V2, V3). This chemistry will be used as decide which barcode whitelist will be used.

    chemistry=\$(
    
    for FASTQ in "\$result"; do
        echo -n ''
        gunzip -cd "\$FASTQ" |
        head -n 400000 | # the first 100k reads should suffice
        awk 'FNR%4==2' | # extract the read sequence
        cut -c1-14 | # the barcode should be in the first 14 bases (we ignore that V2/3 chemistries have longer barcodes)
        awk '
            FNR==1{ fileno++ } # keep track of the file number being processed
            fileno<=3 { version=fileno; barcodes[version","\$0] } # cache barcodes in memory
            fileno==4 { for (version=1; version<=3; version++) if (version","\$0 in barcodes) count[version]++; total++ } # count matches for each chemistry
            END { for (version=1; version<=3; version++) if (count[version]/total > 0.85) printf "V"version } # output chemistries with >85% matches
        ' <(echo "\$V1_BARCODES") <(echo "\$V2_BARCODES") <(echo "\$V3_BARCODES") /dev/stdin
        echo 
    done)

    protocols=\$(tr V v <<<"10x\$chemistry")
    
    echo "\$protocols"
    echo "\$chemistry"
    
    # Based on the chemistry to choose the correct 10x barcide whitelist
    gzip  -dcf \$PWD/whitelist/10x_"\$chemistry"_barcode_whitelist.txt.gz > whitelist.uncompressed.txt
    

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
        -c "\$protocols" \\
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
