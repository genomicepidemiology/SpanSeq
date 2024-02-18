rule list_files:
    input:
        in_dir = "%s" % (config["general"]["sample"]["data_path"])
    output:
        out_file = "%s{sample}%s" % (input_folder, fsa_ext)
    benchmark:
        "benchmarks/{sample}.resumefiles.benchmark.txt"
    run:
        import os
        with open(output.out_file, "w") as out_write:
            for file in os.listdir(input.in_dir):
                file_path = os.path.join(input.in_dir, file)
                if os.path.isfile(file_path):
                    abs_path = os.path.abspath(file_path)
                    out_write.write(abs_path)
                    out_write.write("\n")
