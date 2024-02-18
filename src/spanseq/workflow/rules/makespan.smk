if config["clustering"]["activate"]:
    field_key = 3
elif config["reducing"]["activate"]:
    field_key = 2
else:
    raise ValueError("No clustering has been used before running makespan")

rule ccphylo_makespan:
    input:
        output_dbscan = "{sample}_clusters.tsv"
    output:
        output_makespan = "{sample}_makespan.tsv",
        stats_makespan = "{sample}_makespan_stats.tsv"
    threads:
        config["clustering"]["threads"]
    benchmark:
        "benchmarks/{sample}.makespan.benchmark.txt"
    params:
        machines = config["software"]["makespan"]["machines"],
        method = config["software"]["makespan"]["makespan_method"],
        weight_method = config["software"]["makespan"]["weight_method"],
        field_cluster = field_key,
        ccphylo_path = config["software"]["makespan"]["path"],
        extra_args = config["software"]["makespan"]["extra_args"]
    log:
        "log/{sample}_makespan.log"
    shell:
        "{params.ccphylo_path}ccphylo makespan -i {input.output_dbscan} "\
        "-o {output.output_makespan} -k {params.field_cluster} {params.extra_args}" \
        "-m {params.method} -w {params.weight_method} -l {params.machines} "\
        " > {output.stats_makespan} 2>&1"
