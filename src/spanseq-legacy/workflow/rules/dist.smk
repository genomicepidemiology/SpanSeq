##TODO: fix -m

rule kma_dist:
    input:
        output_index = expand("tmp/{{sample}}{ext}", ext=INDEX_KMA_EXT)
    output:
        output_dist = "tmp/{sample}.phy"
    threads:
        config["clustering"]["threads"]
    benchmark:
        "benchmarks/{sample}.dist.benchmark.txt"
    params:
        memory_disk = config["software"]["kma"]["memory_disk"],
        temp_files = config["general"]["tmpdir"],
        method = config["software"]["kma"]["dist_method"],
        kma_path = config["software"]["kma"]["path"],
        extra_args = config["software"]["kma"]["extra_args_dist"],
        tmpdir = config["general"]["tmpdir"]
    log:
        "log/{sample}_kmadist.log"
    shell:
        "{params.kma_path}kma dist -t_db {params.tmpdir}{sample} -o {output.output_dist} " \
        "-d {params.method} -tmp {params.temp_files} -t {threads} {params.extra_args}"
