import yaml
import os
import sys
import warnings
from .file_mixin import FileValidation
from .pipeline_conf import Pipelines


class Config(dict):

    def __init__(self, configfile=None, envfile=None):

        if configfile is None:
            self.configfile = Config.get_empty_configfile()
        else:
            self.configfile = configfile
        if envfile is None:
            self.envfile = Config.get_empty_envfile()
        else:
            self.envfile = envfile

    @staticmethod
    def get_empty_envfile():
        prev_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        env_file = "{}/config/env_empty.yaml".format(prev_dir)
        env_file = FileValidation.is_valid_file(file=env_file)
        with open(env_file) as file:
            environment = yaml.full_load(file)
            return environment

    @staticmethod
    def get_empty_configfile():
        prev_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        config_file = "{}/config/config_empty.yaml".format(prev_dir)
        conf_file = FileValidation.is_valid_file(file=config_file)
        with open(conf_file) as file:
            configuration = yaml.full_load(file)
            return configuration

    def __setitem__(self, i, y):
        if i == "general":
            self.configfile["general"].update(y)
        elif i == "clustering":
            self.configfile["clustering"].update(y)
        elif i == "software":
            self.configfile["software"].update(y)
        else:
            raise KeyError("""The key {} is invalid. Only keys accepted are"""
                           """ 'general', 'clustering' and 'methods'""".format(
                           i))

    def __getitem__(self, y):
        if i == "general":
            value = self.general_params
        elif i == "clustering":
            value = self.clustering_params
        elif i == "software":
            value = self.clustering_params
        else:
            raise KeyError("""The key {} is invaldi. Only keys accepted are"""
                           """ 'general', 'clustering' and 'methods'""".format(
                           i))
        return value

    def set_params(self, arg_params):

        self.set_general_params(arg_params)
        if arg_params.action == "split":
            self.set_clustering_params(arg_params)
        elif arg_params.action == "reduce":
            self.set_reducing_params(arg_params)

        if arg_params.action is not None:
            scheme = Pipelines.get_pipeline(mode=arg_params.action)
            Pipelines.check_bioseq(scheme=scheme, bioseq=arg_params.seqtype)
        else:
            sys.exit("""In order to run the SpanSeq pipeline, a pipeline """
                     """scheme (-p), config file (-c) or a list of methods """
                     """(-b) has to be selected.""")
        self.set_methods(scheme, arg_params)


    def set_general_params(self, args):
        output_dir = os.path.join(os.path.abspath(args.output_folder), "")
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        self.configfile["general"]["basedir"] = output_dir
        self.configfile["general"]["outfmt"] = args.output_format
        if args.temp_files is None:
            tmp_file = "{}tmp/".format(output_dir)
            if not os.path.exists(tmp_file):
                os.mkdir(tmp_file)
        else:
            tmp_file = args.temp_files
        self.configfile["general"]["resultsdir"] = output_dir
        self.configfile["general"]["tmpdir"] = tmp_file
        self.configfile["general"]["logdir"] = "{}log/".format(output_dir)
        self.configfile["general"]["keep_tmp"] = args.keep_tmp
        self.set_sample_params(args)

    def set_sample_params(self, args):
        self.configfile["general"]["sample"] = dict()
        if (args.input_fasta is None and args.input_folder is None
            and args.input_batch is None):
            raise ValueError("No Input fasta or input folder has been selected.")
        elif (args.input_fasta is not None and args.input_folder is not None
            and args.input_batch is not None):
            raise ValueError("Input fasta and input folder has been selected.")
        else:
            if args.input_fasta is not None:
                input_data = args.input_fasta
                self.configfile["general"]["sample"]["data_format"] = "File"
            elif args.input_batch is not None:
                input_data = args.input_batch
                self.configfile["general"]["sample"]["data_format"] = "Batch"
            else:
                input_data = args.input_folder
                self.configfile["general"]["sample"]["data_format"] = "Folder"
        if not os.path.exists(input_data):
            sys.exit("File {} does not exists".format(args.input_fasta))
        else:
            input_data = os.path.abspath(input_data)
        data_dir = os.path.join(os.path.abspath(os.path.dirname(input_data)),
                                "")
        self.configfile["general"]["datadir"] = data_dir
        filename = os.path.basename(os.path.abspath(input_data))
        self.configfile["general"]["sample"]["data_path"] = os.path.abspath(input_data)
        _name, _fsa_ext = os.path.splitext(filename)
        name = _name.strip()
        fsa_ext = _fsa_ext.strip()
        if fsa_ext != "":
            self.configfile["general"]["sample"]["fsa_ext"] = fsa_ext
            self.configfile["general"]["sample"]["name"] = name
        else:
            self.configfile["general"]["sample"]["fsa_ext"] = '.txt'
            self.configfile["general"]["sample"]["name"] = name

        self.configfile["general"]["sample"]["bioseq_type"] = args.seqtype

    def set_reducing_params(self, args):
        try:
            dist_value = float(args.min_dist)
        except ValueError:
            raise ValueError("""The value of minimum distance between """
                             """samples in different clusters has to be """
                             """in float format, instead of {}""".format(
                             args.min_dist))
        self.configfile["reducing"]["dist_value"] = dist_value

        self.configfile["reducing"]["threads"] = int(args.threads)

    def set_clustering_params(self, args):
        try:
            dist_value = float(args.min_dist)
        except ValueError:
            raise ValueError("""The value of minimum distance between """
                             """samples in different clusters has to be """
                             """in float format, instead of {}""".format(
                             args.min_dist))
        self.configfile["clustering"]["dist_value"] = dist_value

        self.configfile["clustering"]["threads"] = int(args.threads)
        if args.approach != "all" and not args.hobohm1_distance:
            raise ValueError("""When using modes 'hobohm_reduce' and 'hobohm_split', """
                    """a value must be set for -hd/--hobohm1_distance.""")
        else:
            self.configfile["clustering"]["mode"] = args.approach


    def set_software(self, methods, methods_params):
        for method, method_value in methods.items():
            if method_value["run"] == "required":
                self.configfile["software"][method_value["method"]]["activate"] = True
            elif method_value["run"] == "option":
                if method in methods_params:
                    self.configfile["software"][method_value["method"]]["activate"] = True
                    if method_value["param"] is not None:
                        self.configfile["software"][method_value["method"]]["dist_method"] = method_value["param"]
                        self.configfile["clustering"]["dist_value"]/=method_value["factor_dist"]
            else:
                raise KeyError("The method {} is not set well in the pipeline.schemes.yaml".format(method))

    def set_reducer_methods(self, scheme, args):
        self.configfile["reducing"]["activate"] = True
        # Set distance
        for function, methods in scheme["steps"].items():
            if isinstance(methods, dict):
                self.set_software(methods=methods,
                            methods_params=[])
            else:
                raise TypeError(
                    "The Pipeline scheme template is not formatted adequately")
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["kma"],
                                    key="hobohm1_value", value=args.min_dist)


    def set_cluster_methods(self, scheme, args):
        self.configfile["clustering"]["activate"] = True

        # Set distance
        for function, methods in scheme["steps"].items():
            if isinstance(methods, dict):
                self.set_software(methods=methods,
                            methods_params=[args.distanceMethod])
            else:
                raise TypeError(
                    "The Pipeline scheme template is not formatted adequately")

        memory_disk = bool(args.memory_disk)
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["dbscan", "kma"],
                                    key="memory_disk", value=memory_disk)


        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["kma"],
                                    key="hobohm1_value", value=args.hobohm1_distance)


    def set_methods(self, scheme, args):

        if args.action == "split":
            self.set_cluster_methods(scheme, args)
        elif args.action == "reduce":
            self.set_reducer_methods(scheme, args)


        temp_files = str(args.temp_files)
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["dbscan", "kma"],
                                    key="temp_files", value=temp_files)
        kmer_size = int(args.kmer_size)
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["mash", "kma"],
                                    key="kmer_size", value=kmer_size)
        if args.minimizer_size:
            minimizer_size = int(args.minimizer_size)
        else:
            minimizer_size = args.minimizer_size

        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["kma"],
                                    key="minimizer_size", value=minimizer_size)

        prefix = args.prefix
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["kma"],
                                    key="prefix", value=prefix)
        if args.MegaDB:
            megadb = "-ME"
        else:
            megadb = ""
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["kma"],
                                    key="megadb", value=megadb)


        sketch_size = Config.set_sketch_size(configfile=self.configfile,
                                             max_len=args.max_length)
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["mash"],
                                    key="sketch_size", value=sketch_size)
        if isinstance(args.bins, list):
            try:
                machines = [int(x) for x in args.bins]
            except:
                ValueError("""The value in the list for the machines have to"""
                           """ be a list of integers""")
        else:
            try:
                machines = int(args.bins)
            except ValueError:
                raise ValueError("""The value for the machines has to be an"""
                                 """ integer or a list of integers. {} has """
                                 """been used""".format(args.bins))
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["makespan"],
                                    key="machines", value=machines)

        mspan_method = str(args.makespanProcess)
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["makespan"],
                                    key="makespan_method", value=mspan_method)

        w_method = str(args.makespanWeights)
        Config.set_method_subparam(configfile=self.configfile,
                                    method_configs=["makespan"],
                                    key="weight_method", value=w_method)
        

    @staticmethod
    def set_method_subparam(configfile, method_configs, key, value):
        for method in method_configs:
            if bool(configfile["software"][method]["activate"]):
                if value is not None:
                    configfile["software"][method][key] = value
                else:
                    raise ValueError("""The method {method} requires the """
                                     """parameter {key} to be set""".format(
                                     method=method, key=key))

    @staticmethod
    def set_sketch_size(configfile, max_len):
        sketch_size = 1
        if not configfile["software"]["mash"]["activate"]:
            sketch_size = None
        elif max_len is None:
            raise ValueError("""The method Mash requires the """
                             """parameter max length to be set""")
        else:
            while sketch_size < max_len:
                sketch_size+=sketch_size
        return sketch_size


    @staticmethod
    def assign_path(configfile, method, executable, path, conda=True):
        if conda:
            if os.path.exists(path):
                if os.path.isfile(path):
                    dir_path = os.path.dirname(path)
                else:
                    dir_path = path + "/"
                configfile["software"][method]["path"] = dir_path
            else:
                sys.exit("The folder {} does not exists".format(path))
        else:
            if path:
                if os.path.exists(path):
                    if os.path.isfile(path):
                        dir_path = os.path.dirname(path)
                    else:
                        dir_path = path + "/"
                    configfile["software"][method]["path"] = dir_path
                else:
                    sys.exit("The folder {} does not exists".format(dir_path))
            else:
                sp_folder = os.path.dirname(os.path.dirname(
                                                os.path.realpath(__file__)))
                bin_folder = "{}/../../bin/".format(sp_folder)
                exec_file = "{0}/{1}/".format(bin_folder, executable)
                if os.path.exists(exec_file):
                    configfile["software"][method]["path"] = exec_file
                else:
                    executable = Executable(name=executable)
                    repo_dir = executable.install(bin_folder=bin_folder,
                                                  compile=True)
                    configfile["software"][method]["path"] = "{}".format(
                                                                      repo_dir)


    def set_environment(self):
        for method in list(self.envfile["dependencies"].keys()):
            if method not in self.configfile["software"]:
                del self.envfile["dependencies"][method]

    def set_executables(self, args):
        for method in self.configfile["software"]:
            if method == "kma":
                if args.kmaPath:
                    Config.assign_path(configfile=self.configfile,
                                       executable="kma", method=method,
                                       path=os.path.abspath(args.kmaPath))
            elif method == "mash":
                if args.mashPath:
                    Config.assign_path(configfile=self.configfile,
                                       executable="mash", method=method,
                                       path=os.path.abspath(args.mashPath))
            elif method == "dbscan" or method == "makespan":
                if args.ccphyloPath:
                    Config.assign_path(configfile=self.configfile,
                                       method=method, executable="ccphylo",
                                       path=os.path.abspath(args.ccphyloPath))
            else:
                sys.exit("The method {} is not added at SpanSeq".format(method))



    def write_config_file(self, out_file):
        with open(out_file, 'w') as outfile:
            yaml.dump(self.configfile, outfile, default_flow_style=False)
