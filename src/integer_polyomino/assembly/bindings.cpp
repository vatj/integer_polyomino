#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "assembly_api.hpp"

namespace py = pybind11;


PYBIND11_MODULE(assembly, m) {
    m.doc() = R"pbdoc(
        Integer Polyomino Stochastic Assembly Module
        ---------------------------------------------
        .. currentmodule:: assembly
        .. autosummary::
           :toctree: _generate
    )pbdoc"; // optional module docstring

    // m.doc() = R"pbdoc(Integer Polyomino Stochastic Assembly Module)pbdoc";; // optional module docstring

    // Assembly Functions
    m.def("AssemblePlasticGenotype", &AssemblePlasticGenotypeAPI,
      R"pbdoc(A function which builds a genome, create a new phenotype table and return the Phenotype_IDs of the built polyominoes)pbdoc",
      py::arg("genome"), py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 40,
      py::arg("fixed_table") = true, py::arg("determinism") = 1,
      py::arg("table_file") = std::string("None"));

    m.def("AssemblePlasticGenotypeFrequency", &AssemblePlasticGenotypeFrequencyAPI,
      R"pbdoc(A function which builds a genome, create a new phenotype table and return the Phenotype_IDs of the built polyominoes)pbdoc",
      py::arg("genome"), py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 40,
      py::arg("fixed_table") = true, py::arg("determinism") = 1,
      py::arg("table_file") = std::string("None"));

    m.def("AssemblePlasticGenotypes", &AssemblePlasticGenotypesAPI,
      R"pbdoc(A function which builds genomes, create a new phenotype table and return the Phenotype_IDs of the built polyominoes)pbdoc",
      py::arg("genomes"), py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 40,
      py::arg("fixed_table") = true, py::arg("determinism") = 1,
      py::arg("table_file") = std::string("None"));

    m.def("AssemblePlasticGenotypesFrequency", &AssemblePlasticGenotypesFrequencyAPI,
      R"pbdoc(A function which builds genomes, create a new phenotype table and return the Phenotype_IDs of the built polyominoes and their frequencies)pbdoc",
      py::arg("genomes"), py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 40,
      py::arg("fixed_table") = true, py::arg("determinism") = 1,
      py::arg("table_file") = std::string("None"));

    #ifdef VERSION_INFO
      m.attr("__version__") = VERSION_INFO;
    #else
      m.attr("__version__") = "dev";
    #endif
}
