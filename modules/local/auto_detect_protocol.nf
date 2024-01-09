process AUTO_DETECT_PROTOCOL {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::jq=1.6'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.6' :
        'biocontainers/jq:1.6' }"

    input:
    // the first FastQ file in `reads` is expected to contain the cell barcodes
    tuple val(meta), path(reads)
    val aligner
    path projectDir

    output:
    tuple val(meta), path(reads), emit: ch_fastq
    env PROTOCOL, emit: protocol
    env EXTRA_ARGS, emit: extra_args
    path "*.gz", emit: whitelist, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    # convert protocols.json to table
    TABLE=\$(
        jq -r '
            ."$aligner" |
            to_entries[] |
            "\\(.key)\\t\\(.value.protocol//"")\\t\\(.value.whitelist//"")\\t\\(.value.extra_args//"")"
        ' "${projectDir}/assets/protocols.json"
    )

    # iterate over all protocols defined for the selected aligner
    KEY=\$(cut -f1 <<<"\$TABLE" | while read KEY; do

        # uncompress whitelist
        WHITELIST=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f3)
        [ -n "\$WHITELIST" ] || continue # skip protocols without whitelist
        gzip -dcf "${projectDir}/\$WHITELIST" > barcodes

        # subsample the FastQ file
        gzip -dcf "${reads[0]}" |
        awk 'FNR % 4 == 2' | # extract the read sequence from FastQ
        head -n 100000 > reads || true # the first 100k reads should suffice

        # extract the barcodes from the FastQ reads and count how many are valid barcodes
        awk -v KEY="\$KEY" '
            { \$0 = substr(\$0, 1, 14) } # the barcode is in the first 14 bases; 10X V2/3 barcodes are trimmed
            FILENAME == "barcodes" { barcodes[\$0] } # cache barcodes in memory
            FILENAME == "reads" && \$0 in barcodes { count++ } # count matches for each chemistry
            END { if (count/FNR > 0.85) print KEY } # output chemistries with >85% matches
        ' barcodes reads

    done)

    if [ \$(wc -w <<<"\$KEY") -ne 1 ]; then
         echo "ERROR: protocol detection failed: \$KEY"
         exit 1
    fi

    # extract attributes of chosen protocol
    PROTOCOL=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f2)
    WHITELIST=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f3)
    [ -z "\$WHITELIST" ] || cp "${projectDir}/\$WHITELIST" .
    EXTRA_ARGS=\$(grep -w "^\$KEY" <<<"\$TABLE" | cut -f4)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jq: \$(jq --version | cut -d- -f2)
    END_VERSIONS
    """
}
