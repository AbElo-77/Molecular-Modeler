"""
Setup script for building the molecular_interactions C++ extension module
"""
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_dir = os.path.dirname(os.path.abspath(__file__))

ext_modules = [
    Pybind11Extension(
        "molecular_interactions",
        [
            "pybind11_module.cpp",
            "radii_lists.cpp",
            "nonbond_interactions.cpp",
            "bonded_interactions.cpp",
        ],
        include_dirs=[ext_dir],
        extra_compile_args=['-O3'] if sys.platform != 'win32' else ['/O2'],
        language='c++'
    ),
]

setup(
    name="molecular_interactions",
    version="1.0.0",
    author="Molecular Modeler Team",
    description="High-performance C++ molecular interaction calculations",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.8",
)
