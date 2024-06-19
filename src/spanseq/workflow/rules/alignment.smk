input_folder = datadir
if not config["general"]["sample"]["data_format"] == "File":
    raise TypeError("Only 'File' data format is available for alignment at the moment")

rule alignment_dist:
    input:
        in_fsa = "%s{sample}%s" % (input_folder, fsa_ext)
    output:
        dist = "tmp/{sample}.phy"
    log:
        "log/{sample}_aln.log"
    params:
        bioseq = config["general"]["sample"]["bioseq_type"],
        data_dir = config["general"]["datadir"],
        max_len = config["software"]["ggsearch36"]["max_length"],
        ggsearch_path = config["software"]["ggsearch36"]["path"],
        extra_args_aln = config["software"]["ggsearch36"]["extra_args_aln"],
        script_folder = SCRIPTFOLDER
    threads:
        config["clustering"]["threads"]
    benchmark:
        "benchmarks/{sample}.aln.benchmark.txt"
    shell:
        """
        python {params.script_folder}/wrapper_aln.py -I -p {threads} -m {params.max_len} -s {params.bioseq} -a ggsearch36 -o tmp/ -t {tmpdir}/aln_files -f {input.in_fsa}
        """
