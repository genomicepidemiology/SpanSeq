rule merged_df:
    input: 
        clstr = "{sample}_clusters.tsv",
        makespan = "{sample}_makespan.tsv",
    output:
        df_merged = "{sample}_partitions.tsv"
    log:
        "log/{sample}_merged_df.log"
    benchmark:
        "benchmarks/{sample}.merged_df.benchmark.txt"
    params:
        script = "%s/merge_df.py" % SCRIPTFOLDER,
        outdir = config["general"]["resultsdir"]
    shell:
        "python {params.script} -c {input.clstr} -m {input.makespan} -o {output.df_merged}"

rule create_fastas:
    input:
        df_merged = "{sample}_partitions.tsv",
        fasta = "%s{sample}.fsa" % datadir
    output:
        fastas_cluster = expand("{{sample}}_M{bins}.fsa",
                                bins=LIST_FSA),
    log:
        "log/{sample}_create_fasta.log"
    benchmark:
        "benchmarks/{sample}.create_fasta.benchmark.txt"
    params:
        script = "%s/join_clusters.py" % SCRIPTFOLDER,
        bins = config["clustering"]["machines"],
        outdir = config["general"]["resultsdir"]
    shell:
        "python {params.script} -d {input.df_merged} -o {params.outdir} -f {input.fasta} -b {params.bins}"
