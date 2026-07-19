"""Build the authoritative ddca module as a Cython extension."""

from Cython.Build import cythonize
from setuptools import setup


setup(
    ext_modules=cythonize(
        ["ddca.py"],
        compiler_directives={"language_level": "3"},
    )
)
