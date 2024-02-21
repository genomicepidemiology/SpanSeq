from os import path
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

from src.spanseq import __version__, __author__, __email__, __copyright__
from src.spanseq.utils import manage_executables, environment_mixin


this_directory = path.abspath(path.dirname(__file__))


def friendly(command_subclass):
    """A decorator for classes subclassing one of the setuptools commands.

    It modifies the run() method so that it prints a friendly greeting.
    """
    orig_run = command_subclass.run
    req_file = "{}/data/executable_requirements.yml".format(this_directory)
    bin_dir = "{}/bin/".format(this_directory)

    def modified_run(self):
        environment_mixin.Environment.get_envname()
        orig_run(self)
        installing = manage_executables.Build_Ext(requirements_file=req_file,
                                                  bin_folder=bin_dir)
        installing.check_executables()

    command_subclass.run = modified_run
    return command_subclass

...

@friendly
class CustomDevelopCommand(develop):
    pass

@friendly
class CustomInstallCommand(install):
    pass

#class PostDevelopCommand(develop):
#    """Pre-installation for development mode."""
#    def run(self):
#        #check_call("apt-get install this-package".split())
#        develop.run(self)

#class PostInstallCommand(install):
#    """Pre-installation for installation mode."""
#    def run(self):
        #check_call("apt-get install this-package".split())
#        install.run(self)


#class PostCommand(install,develop):
#    """Post-installation for installation mode."""
#    def run(self):
#        install.run(self)

# read the contents of your README file


with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="spanseq",
    version=__version__,
    cmdclass={
        'develop': CustomDevelopCommand,
        'install': CustomInstallCommand,},
    #cmdclass=versioneer.get_cmdclass(),
    url="https://bitbucket.org/genomicepidemiology/spanseq.git",
    author=__author__,
    author_email=__email__,
    zip_safe=False,
    description="SpanSeq - workflows for homology division of databases of biological sequences.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires='>=3.8, <3.11',
    package_data={
        "": [
            "src/spanseq/*",
        ]
    },
    data_files=[(".", ["README.md", "LICENSE.txt"])],
    include_package_data=True,
    install_requires=["pandas>=1.2, <1.5", "biopython"],
    entry_points={"console_scripts": ["spanseq = spanseq.spanseq:main"]},
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
    packages=find_packages(include=['spanseq', 'spanseq.utils']),
    package_dir={'':'src'}
)
