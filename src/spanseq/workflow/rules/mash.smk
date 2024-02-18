input_folder = datadir
if config["general"]["sample"]["data_format"] == "File":
    param_input = "-i"
else:
    param_input = "-l"
    if config["general"]["sample"]["data_format"] == "Folder":
        input_folder = tmpdir
        include: "resume_folder.smk"
        fsa_ext = ".txt"

rule mash_dist:
    input:
        in_fsa = "%s{sample}%s" % (input_folder, fsa_ext)
    output:
        dist = "tmp/{sample}.phy"
    log:
        "log/{sample}_mash.log"
    params:
        bioseq = config["general"]["sample"]["bioseq_type"],
        data_dir = config["general"]["datadir"],
        mash_path = config["software"]["mash"]["path"],
        sketch_size = config["software"]["mash"]["sketch_size"],
        kmer_size = config["software"]["mash"]["kmer_size"],
        extra_args_sketch = config["software"]["mash"]["extra_args_sketch"],
        extra_args_triangle = config["software"]["mash"]["extra_args_triangle"]
    threads:
        config["clustering"]["threads"]
    benchmark:
        "benchmarks/{sample}.mash.benchmark.txt"
    shell:
        """
        bioseq=("{params.bioseq}")
        if [ "$bioseq" = "nucleotides" ]
        then
            {params.mash_path}mash triangle {param_input} {input.in_fsa} \
            -k {params.kmer_size} -s {params.sketch_size} -i  \
            -p {threads} {params.extra_args_sketch} > {output.dist} \
            {params.extra_args_triangle}
        else
            {params.mash_path}mash triangle {param_input} {input.in_fsa} \
            -k {params.kmer_size} -s {params.sketch_size} -i -a \
            -p {threads} {params.extra_args_sketch} > {output.dist}\
             {params.extra_args_triangle}
        fi
        """
