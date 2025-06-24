from git import Repo
import subprocess
import os
import sys

class CondaError(Exception):
    """Exception raised for errors in the input salary.

    Attributes:
        salary -- input salary which caused the error
        message -- explanation of the error
    """

    def __init__(self, env, message=("The conda environment {} must be "
                                     "activated for this step")):
        self.message = message.format(env)
        super().__init__(self.message)

class Executable:

    repositories = {
        "ccphylo": {
            "repo":"https://bitbucket.org/genomicepidemiology/ccphylo.git",
            "branch": "makespan"}
        }

    def __init__(self, name, webpath=None):
        self.name = name
        self.platform_install = None
        if name in Executable.repositories:
            self.platform_install = "git"
            self.address = Executable.repositories[name]["repo"]
        elif webpath is not None:
            self.address = webpath
        else:
            raise OSError("""The executable {} needs a webpath or repository"""
                          """""".format(name))

    def install(self, bin_folder, compile=False):
        repo_dir = "{}/{}".format(bin_folder, self.name)
        Repo.clone_from(self.address, repo_dir,
                        branch=Executable.repositories[self.name]["branch"])
        if compile:
            subprocess.call("make", shell=True, cwd=repo_dir)
        return repo_dir


class Environment:

    standard_names = ["spanseq-light", "spanseq"]

    def get_envname():

        try:
            conda_env = os.environ["CONDA_DEFAULT_ENV"]
        except KeyError:
            raise CondaError(env="spanseq")
            #raise CondaError(("ENVIRONMENT ERROR: No conda environment has been "
            #    "activated. SpanSeq runs with the environments from "
            #    "spanseqenv_light.yml (spanseq-light) or spanseq.yml "
            #    "(spanseq)."))
        if os.path.split(conda_env)[1] not in Environment.standard_names:
            sys.exit(("SpanSeq requires that the conda environment activated"
                " is spanseq-light (spanseqenv_light.yml) or spanseq "
                "(spanseq.yml) not {}".format(os.path.split(conda_env)[1])))
        else:
            return conda_env

    def check_version():
        subprocess.check_output(["conda", "list"])
