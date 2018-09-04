import os
import sys
from setuptools import setup, PEP420PackageFinder, Extension
from distutils.command.build_ext import build_ext
from distutils.ccompiler import new_compiler, get_default_compiler
import sysconfig


class LocalBuild(build_ext):
    """A builder that knows our c++ compiler."""

    def initialize_options(self):
        super().initialize_options()
        # self.compiler = new_compiler(plat=sys.platform, verbose=1)
        # print(f"my compiler {self.compiler}")

    def finalize_options(self):
        super().finalize_options()
        print(f"compiler {type(self.compiler)}")

    def run(self):
        super().run()
        for n, v in self.compiler.__dict__.items():
            print(f"compiler {n}: {v}")


extra_compile_args = sysconfig.get_config_var("CFLAGS").split()
extra_compile_args.append("--std=c++11")
lifetable = Extension(
    "fundem._lifetable",
    sources=["src/lifetable.cpp"],
    include_dirs=["include"],
    libraries=["gsl"],
    language="c++",
    extra_compile_args=extra_compile_args,
)

setup(
    name="fundem",
    version="0.0.1",
    packages=PEP420PackageFinder.find("src"),
    package_data={},
    package_dir={"": "src"},
    install_requires=["numpy", "scipy", "toml"],
    extras_require={
        "testing": ["hypothesis", "pytest", "pytest-mock"],
        "documentation": ["sphinx", "sphinx_rtd_theme", "sphinx-autobuild",
                          "sphinxcontrib-napoleon"],
        "ihme_databases": ["db_tools", "db_queries"],
    },
    ext_modules=[lifetable],
    cmdclass={"build_ext": LocalBuild},
    entry_points={
        "console_scripts": []
    },
    scripts=[],
    zip_safe=False,
    classifiers=[
        "Intendend Audience :: Science/Research",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Statistics",
    ],
)
