import os
from setuptools import setup, find_packages

setup(
    name="abg",
    version="0.0.01",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    packages=["abg"],
    py_modules=[
        "abg/base",
        "abg/bio",
        "abg/compute",
        "abg/matvec, abg/pdb_parser",
        "abg/settings",
    ],
    entry_points={"console_scripts": ["calc_abg = abg.compute:main"]},
    include_package_data=True,
)
