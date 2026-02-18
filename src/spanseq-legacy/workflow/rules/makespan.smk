if config["clustering"]["activate"]:
    field_key = 3
elif config["reducing"]["activate"]:
    field_key = 2
else:
    raise ValueError("No clustering has been used before running makespan")

if config["software"]["makespan"]["imbalance_file"]:
    input_makespan = "{sample}_clusters_class.tsv"
    class_imbalance = True
else:
    input_makespan = "{sample}_clusters.tsv"
    class_imbalance = False

rule add_classes:
    input:
        output_dbscan = "{sample}_clusters.tsv"
    output:
        input_makespan = "{sample}_clusters_class.tsv"
    benchmark:
        "benchmarks/{sample}.add_class.benchmark.txt"
    params:
        file_class = config["software"]["makespan"]["imbalance_file"],
        script = "%s/add_class.py" % SCRIPTFOLDER,
    log:
        "log/{sample}_add_class.log"
    shell:
        "python {params.script} -d {input.output_dbscan} -o {output.input_makespan} -c {params.file_class}"

rule ccphylo_makespan:
    input:
        output_dbscan = input_makespan
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
        extra_args = config["software"]["makespan"]["extra_args"],
        classes = config["software"]["makespan"]["class_columns"],
    log:
        "log/{sample}_makespan.log"
    shell:
        """
        imb_classes=("{params.classes}")
        if [ "$imb_classes" = False ]
        then
            {params.ccphylo_path}ccphylo makespan -i {input.output_dbscan} -o {output.output_makespan} -k {params.field_cluster} {params.extra_args} -m {params.method} -w {params.weight_method} -l {params.machines} > {output.stats_makespan} 2>&1
        else
            {params.ccphylo_path}ccphylo makespan -i {input.output_dbscan} -o {output.output_makespan} -k {params.field_cluster} {params.extra_args} -m {params.method} -w {params.weight_method} -l {params.machines} -c {params.classes} > {output.stats_makespan} 2>&1
        fi
        """