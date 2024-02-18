## TODO: Add -H option
rule ccphylo_dbscan:
    input:
        output_dist = "tmp/{sample}.phy"
    output:
        output_dbscan = "{sample}_clusters.tsv"
    threads:
        config["clustering"]["threads"]
    benchmark:
        "benchmarks/{sample}.dbscan.benchmark.txt"
    params:
        dist_value = config["clustering"]["dist_value"],
        memory_disk = config["software"]["dbscan"]["memory_disk"],
        temp_files = config["general"]["tmpdir"],
        ccphylo_path = config["software"]["dbscan"]["path"]
    log:
        "log/{sample}_dbscan.log"
    shell:
        "{params.ccphylo_path}ccphylo dbscan -i {input.output_dist} " \
        "-e {params.dist_value}  -p "\
        "-T {params.temp_files} -o {output.output_dbscan}"
