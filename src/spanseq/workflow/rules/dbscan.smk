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
        mem_flag = "-H" if config["software"]["dbscan"]["memory_disk"] else "",
        temp_files = config["general"]["tmpdir"],
        ccphylo_path = config["software"]["dbscan"]["path"]
    log:
        "log/{sample}_dbscan.log"
    shell:
        "{params.ccphylo_path}ccphylo dbscan -i {input.output_dist} " \
        "-e {params.dist_value} {params.mem_flag} -p " \
        "-T {params.temp_files} -o {output.output_dbscan}"
