import os
import yaml
from .file_mixin import FileValidation


class Pipelines:


    @staticmethod
    def get_pipelines_schemes():
        prev_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        empty_config = "{}/config/pipelines.schemes.yaml".format(prev_dir)
        return empty_config

    @staticmethod
    def load_pipeline_config(pipeline_config):
        pipe_conf_file = FileValidation.is_valid_file(file=pipeline_config)
        with open(pipe_conf_file) as file:
            pipelines = yaml.full_load(file)
            return pipelines

    @staticmethod
    def get_pipeline(mode):
        config_file = Pipelines.get_pipelines_schemes()
        pipeline_config = Pipelines.load_pipeline_config(config_file)
        try:
            scheme = pipeline_config["pipelines"][str(mode)]
        except KeyError:
            raise KeyError("The scheme {} does not exists".format(mode))
        return scheme

    @staticmethod
    def check_bioseq(scheme, bioseq):
        if bioseq not in scheme["bio_seq"]:
            possible_schemes = []
            pipeline_config = Pipelines.load_pipeline_config(config_file)
            for sch in pipeline_config["schemes"]:
                if bioseq in pipeline_config["schemes"][sch]["bioseq_type"]:
                    possible_schemes.append(
                                    pipeline_config["schemes"][sch]["number"])
            raise TypeError("""The scheme {} cannot be used with sequences """
                            """made of {}. Only the schemes {} can be """
                            """selected.""".format(scheme["number"],
                                                  ", ".join(possible_schemes)))

