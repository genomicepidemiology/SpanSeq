##TODO: fix -m

if config["general"]["sample"]["bioseq_type"] == "nucleotides":
    if config["general"]["reduce"]>=0.9:
        word_size = 8
    elif config["general"]["reduce"]>=0.88 and config["general"]["reduce"]<0.9:
        word_size = 7
    elif config["general"]["reduce"]>=0.85 and config["general"]["reduce"]<0.88:
        word_size = 6
    elif config["general"]["reduce"]>=0.80 and config["general"]["reduce"]<0.85:
        word_size = 5
    elif config["general"]["reduce"]>=0.75 and config["general"]["reduce"]<0.8:
        word_size = 4
else:
    if config["general"]["reduce"]>=0.7:
        word_size = 5
    elif config["general"]["reduce"]>=0.6 and config["general"]["reduce"]<0.7:
        word_size = 4
    elif config["general"]["reduce"]>=0.5 and config["general"]["reduce"]<0.6:
        word_size = 3
    elif config["general"]["reduce"]>=0.4 and config["general"]["reduce"]<0.5:
        word_size = 2


rule kma_dist:
    input:
        in_fsa = "%s{sample}%s" % (input_folder, fsa_ext)
    output:
        reduced_fsa = "tmp/{sample}_cdhit.fsa"
        reduced_cluster = "tmp/{sample}_cdhit.fsa.cltsr"
    threads:
        config["clustering"]["threads"]
    benchmark:
        "benchmarks/{sample}.cdhit.benchmark.txt"
    params:
        data_type = config["general"]["sample"]["bioseq_type"],
        distance_reduce = config["general"]["reduce"],
        tmpdir = config["general"]["tmpdir"],
        wordsize = word_size
    log:
        "log/{sample}_cdhit.log"
    shell:
        """
        data_type=("{params.data_type}")
        if [ "$data_type" = "nucleotides" ]
        then
            cd-hit -i {input.in_fsa} -o {output.reduced_fsa} -n {params.wordsize} -c {params.distance_reduce}
        else
            cd-hit-est -i {input.in_fsa} -o {output.reduced_fsa} -n {params.wordsize} -c {params.distance_reduce}
        fi
        """
