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
        "bin/spdw-bpp-analysis-profile.py",
        "bin/spdw-bpp-create-analyses.py",
        "bin/spdw-bpp-explore-priors.py",
        "bin/spdw-bpp-extract-a10-tree.py",
        "bin/spdw-bpp-fix-traces.py",
        "bin/spdw-calc-inv-gamma-prior.py",
        "bin/spdw-delineate-create-analyses-from-trees.py",
        "bin/spdw-delineate-create-analyses.py",
        "bin/spdw-delineate-evaluate-analysis.py",
        "bin/spdw-fsc26-create-run.py",
        "bin/spdw-gen-seqs.py",
        "bin/spdw-plotcoaltimes.py",
        "bin/spdw-plot-invgamma.py",
        "bin/spdw-sim-bdcoal-trees.py",
        "bin/spdw-sim-bd-trees.py",
        "bin/spdw-sim-coal-trees.py",
        "bin/spdw-sim-protractedspeciation-trees.py",
    ],
    test_suite = "tests",
    url="http://github.com/jeetsukumaran/spdw",
    license="LICENSE.txt",
    description="Species Delimitation Workshop",
    long_description=_read(["README.txt"]),
    # install_requires=[ ],
)
