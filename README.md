# Genotype-Phenotype Map Integer Polyomino

Code repository for the analysis of the Genotype-Phenotype Map Properties of the Integer Polyomino Model.

## Overview
Custom GP-map analysis tools for the Integer Polyomino Model are coded in C++. These tools rely on a core algorithm to build polyominoes from lists of integer hosted in [this repository](https://github.com/ASLeonard/polyomino_core). 
The python API uses [pybind11](https://pybind11.readthedocs.io/en/stable/) to expose the underlying C++ to the python interpreters. A minimal documentation of the python API is provided using Sphinx and hosted on the Theory of Condensed Matter Group servers.
Additionally, example jupyter notebooks are made available to get started with the python API.

## Installation

This module can be installed as any python package :
```shell
pip install git+https://github.com/vatj/integer_polyomino.git
```
Note that compiling requires a relatively recent C++ compiler such as clang++6 or g++7 as well as libraries such as openml

## Quickstart

Two main quickstart jupyter notebooks are provided to examplify the basic functionality of the repository. Other examples notebooks are also present which can be used for more advance use cases.

### Assembly

This notebook contains a minimal example of how to use the polyomino assembly API. It is intended as a tool to map a list of integer to a polyomino, i.e. a genotype to its corresponding phenotype. In practice the function only returns an ID which corresponds to an entry in a phenotype_table file. Each entry is a matrix of boolean value. A simple plot function is provided to map this matrix onto a polyomino shape.

### GP-Map

This notebook contains an example of how to generate the the analysis of a integer polyomino genotype-phenotype map. All genotype contains 2 genes and at most 7 different labels (1 neutral and 3 interacting pairs). The output is written to 2 different files, one for genome-level metrics and one for phenotype-level metrics. These files can be loaded in panda dataframes for easy analysis and plotting.

## Directory Layout

+ docs/ (integration with sphinx)
+ libs/ (git submodule dependencies)
+ notebooks/ (working examples)
+ src/
  + integer_polyomino/
    + assembly/ (Mostly pybind API for custom polyomino_core features)
    + gpmap/ (C++ implementation of the gpmap analysis algorithms + pybind API definition)
    + scripts/ (python utility file such as plotting functions)
+ tests/
  

