#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
from setuptools import setup, find_packages

def _read(path_components, **kwargs):
    path = os.path.join(os.path.dirname(__file__), *path_components)
    if sys.version_info.major < 3:
        return open(path, "rU").read()
    else:
        with open(path, encoding=kwargs.get("encoding", "utf8")) as src:
            s = src.read()
        return s

project_init = _read(["src", "spdw", "__init__.py"])
__version__ = re.match(r".*^__version__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M).group(1)
__project__ = re.match(r".*^__project__\s*=\s*['\"](.*?)['\"]\s*$.*", project_init, re.S | re.M).group(1)

setup(
    name=__project__,
    version=__version__,
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=find_packages("src"),
    package_dir={"": "src"},
    scripts=[
        "bin/calc-inv-gamma-prior.py",
        "bin/compare-coalescent-vs-msc.sh",
        "bin/spdw-build-bpp-jobs.py",
        "bin/spdw-build-delineate-jobs.py",
        "bin/spdw-build-fsc26-run.py",
        "bin/spdw-evaluate-delineate-jobs.py",
        "bin/spdw-extract-bpp-a10-tree.py",
        "bin/spdw-fix-bpp-traces.py",
        "bin/spdw-gen-seqs.py",
        "bin/spdw-plotcoaltimes.py",
        "bin/spdw-plot-invgamma.py",
        "bin/spdw-sim-bdstruct-coaltrees.py",
        "bin/spdw-sim-bdtrees.py",
        "bin/spdw-sim-coaltrees.py",
        "bin/spdw-sim-protractedspeciationtrees.py",
    ],
    test_suite = "tests",
    url="http://github.com/jeetsukumaran/spdw",
    license="LICENSE.txt",
    description="Species Delimitation Workshop",
    long_description=_read(["README.txt"]),
    # install_requires=[ ],
)
