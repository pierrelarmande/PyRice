from setuptools import setup

with open('README.md','r') as fh:
    long_description = fh.read()

setup(
    name = 'pyrice',
    version = '0.0.1',
    description = 'PyRice: a Python package for functional analysis of rice genes',
    py_modules =['multi_query','utils'],
    license = 'MIT',
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    packages=["pyrice"],
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/pierrelarmande/PyRice/tree/dev",
    author = "Quan Do and Pierre Larmande",
    author_email = "dohongquan1612@gmail.com",
    include_package_data = True
)
