import os
import argparse

class FileValidation:

    @staticmethod
    def is_valid_file(file, parser=None):
        if not os.path.exists(file):
            if parser is None:
                raise OSError("The file %s does not exist!" % file)
            else:
                parser.error("The file %s does not exist!" % file)
        else:
            return file

    @staticmethod
    def readable_dir(prospective_dir):
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            return prospective_dir
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))
