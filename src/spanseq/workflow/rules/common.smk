from snakemake.utils import validate
import os

validate(config, schema="../../config/config.schema.yaml")
current_folder = os.path.dirname(os.path.realpath(__file__))

INDEX_KMA_EXT = [".comp.b", ".length.b",".name", ".seq.b"]
LIST_FSA = list(range(1, int(config["software"]["makespan"]["machines"])+1))


my_basedir = workflow.current_basedir
CONDAENV = "{}/../envs".format(my_basedir)
SCRIPTFOLDER = "{}/../scripts".format(my_basedir)
#conda_envfile = config["general"]["env_file"]

sample = config["general"]["sample"]["name"]
fsa_ext = config["general"]["sample"]["fsa_ext"]
datadir = config["general"]["datadir"]
resultsdir = config["general"]["resultsdir"]
tmpdir = config["general"]["tmpdir"]
logdir = config["general"]["logdir"]

#if config["general"]["sample"]["data_format"] == "File":
#    expand_dir = ""
#else:
#    expand_dir = "/*"
