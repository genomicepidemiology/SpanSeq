input_folder = datadir
if config["general"]["sample"]["data_format"] == "File":
    param_input = "-i"
else:
    param_input = "-batch"
    if config["general"]["sample"]["data_format"] == "Folder":
        input_folder = tmpdir
        include: "resume_folder.smk"
        fsa_ext = ".txt"
if config["reducing"]["activate"]:
    output_hobohm = "{sample}_clusters.tsv"
else:
    output_hobohm = "{sample}_hobohmclusters.tsv"


if config["clustering"]["mode"] == "hobohm_reduce" or config["reducing"]["activate"]:
    rule kma_index_hobohm1:
        input:
            in_fsa = "%s{sample}%s" % (input_folder, fsa_ext)
        output:
            output_hobohm = output_hobohm,
            output_index = expand("tmp/{{sample}}{ext}",
                                  ext=INDEX_KMA_EXT),
        threads:
            config["clustering"]["threads"]
        benchmark:
            "benchmarks/{sample}.kmaindex.benchmark.txt"
        params:
            hobhm_value = config["software"]["kma"]["hobohm1_value"],
            kma_path = config["software"]["kma"]["path"],
            extra_args = config["software"]["kma"]["extra_args_index"],
            kmer_size = config["software"]["kma"]["kmer_size"],
            minimizer_size = config["software"]["kma"]["minimizer_size"],
            prefix = config["software"]["kma"]["prefix"]
        log:
            kma = "%s{sample}_kmaindex.log" % logdir
        shell:
            """
            minimizer=("{params.minimizer_size}")
            if [ "$minimizer" = False ]
            then
                {params.kma_path}kma index {param_input} {input.in_fsa} -o {tmpdir}{sample} -k {params.kmer_size} -hq {params.hobhm_value} -ht {params.hobhm_value} -and {params.extra_args} -Sparse - > {output.output_hobohm} 2> {log.kma}
            else
                {params.kma_path}kma index {param_input} {input.in_fsa} -o {tmpdir}{sample} -k {params.kmer_size} -hq {params.hobhm_value} -ht {params.hobhm_value} -and {params.extra_args} -m {params.minimizer_size} -Sparse - > {output.output_hobohm} 2> {log.kma}
            fi
            """
else:
    rule kma_index:
        input:
            in_fsa = "%s{sample}%s" % (input_folder, fsa_ext)
        output:
            output_index = expand("tmp/{{sample}}{ext}",
                                  ext=INDEX_KMA_EXT),
        threads:
            config["clustering"]["threads"]
        benchmark:
            "benchmarks/{sample}.kmaindex.benchmark.txt"
        params:
            kma_path = config["software"]["kma"]["path"],
            extra_args = config["software"]["kma"]["extra_args_index"],
            kmer_size = config["software"]["kma"]["kmer_size"],
            fsa_ext = config["general"]["sample"]["fsa_ext"],
            minimizer_size = config["software"]["kma"]["minimizer_size"],
            prefix = config["software"]["kma"]["prefix"],
            megadb = config["software"]["kma"]["megadb"]
        log:
            kma = "%s{sample}_kmaindex.log" % logdir
        shell:
            """
            minimizer=("{params.minimizer_size}")
            if [ "$minimizer" = False ]
            then
                {params.kma_path}kma index {param_input} {input.in_fsa} -Sparse {params.prefix} {params.megadb} -k {params.kmer_size} -o {tmpdir}{sample} {params.extra_args} 2> {log.kma}
            else
                {params.kma_path}kma index {param_input} {input.in_fsa} -Sparse {params.prefix} {params.megadb} -k {params.kmer_size} -m {params.minimizer_size} -o {tmpdir}{sample} {params.extra_args} 2> {log.kma}
            fi
            """
