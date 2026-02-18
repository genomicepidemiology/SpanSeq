import subprocess
import re
import yaml
import os
import operator


class Build_Ext:
    operators_map = {"==": operator.eq, "!=": operator.ne,
                     ">=": operator.ge, ">": operator.gt,
                     "<=": operator.le, "<": operator.lt}

    def __init__(self, requirements_file, bin_folder):
        with open(requirements_file, 'r') as stream:
            self.requirements_data = yaml.safe_load(stream)
        if os.path.isdir(bin_folder):
            self.bin_fold = bin_folder
        else:
            raise OSError("The bin folder ({}) do not exists".format(
                                                                bin_folder))

    @staticmethod
    def get_version_std(software, command):
        try:
            result = subprocess.run([str(software), str(command)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True)
            return result.stdout

        except subprocess.CalledProcessError as e:
            return None

    @staticmethod
    def get_version(string_version, software):
        match_version1 = re.findall(r"[0-9]+[.][0-9]+[.][0-9]+", string_version)
        match_version2 = re.findall(r"[0-9]+[.][0-9]+", string_version)
        matched_version = None
        for match_version in [match_version1, match_version2]:
            if len(match_version) == 1:
                matched_version = match_version[0]
                break
            elif len(match_version) > 1:
                raise ValueError("More than one version {} for the software {}".format(
                                    ",".join(match_version), software))
        if matched_version is None:
            raise ValueError("No version has been able to be found for the software {} ({})".format(
                                software, string_version))
        else:
            return matched_version

    @staticmethod
    def check_version(software, orig_operator, req_version, current_version):
        for o_n, c_n in zip(req_version.split("."),
                            current_version.split(".")):
            if Build_Ext.operators_map[orig_operator](o_n, c_n):
                continue
            else:
                raise ValueError("The software {} installed".format(software))

    def check_executables(self):
        for software, val in self.requirements_data.items():
            std_version = Build_Ext.get_version_std(software=software,
                                                    command=val["command"])
            version_installed = Build_Ext.get_version(
                                                    string_version=std_version,
                                                    software=software)
