#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "gpmap_api.hpp"

namespace py = pybind11;


PYBIND11_MODULE(gpmap, m) {
    m.doc() = R"pbdoc(
        Integer Polyomino GPMap Analysis Module
        --------------------------------------------
        .. currentmodule:: gpmap
        .. autosummary::
           :toctree: _generate
    )pbdoc"; // optional module docstring

    // Genome Space generator functions

    m.def("MinimalGenomes", &MinimalGenomesAPI,
      "Create a new phenotype table and return generate minimal genomes for (n_genes, colours)",
      py::arg("n_genes") = 2, py::arg("low_colour") = 0, py::arg("high_colour") = 6,
      py::arg("genome_file") = std::string("None"), py::arg("phenotype_builds") = 40,
      py::arg("threshold") = 0.25, py::arg("fixed_table") = true, py::arg("determinism") = 1);

    m.def("MinimalGenomesVoid", &MinimalGenomesVoidAPI,
      "Create a new phenotype table, generate minimal genomes for (n_genes, colours) and save them in genome_file",
      py::arg("n_genes") = 2, py::arg("low_colour") = 0, py::arg("high_colour") = 6,
      py::arg("genome_file") = std::string("None"), py::arg("phenotype_builds") = 40,
      py::arg("threshold") = 0.25, py::arg("fixed_table") = true, py::arg("determinism") = 1);

    // GPmap functions

    m.def("MinimalMap", &MinimalMapAPI,
      "Return the GP map built from a list of genomes or a genome file.",
      py::arg("genomes") = std::vector <Genotype> (), py::arg("table_file") = std::string("None"),
      py::arg("gpmap_file") = std::string("None"), py::arg("genome_file") = std::string("None"),
      py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 250,
      py::arg("fixed_table") = true, py::arg("determinism") = 1);

    m.def("MinimalMapVoid", &MinimalMapVoidAPI,
      "Save the GP map built from a list of genomes or a genome file.",
      py::arg("genomes") = std::vector <Genotype> (), py::arg("table_file") = std::string("None"),
      py::arg("gpmap_file") = std::string("None"), py::arg("genome_file") = std::string("None"),
      py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 250,
      py::arg("fixed_table") = true, py::arg("determinism") = 1);

    // Neighbourhood

    m.def("GenotypeNeighbourhood", &GenotypeNeighbourhoodAPI,
      "Return the list of the single mutation neighbour of a genome",
      py::arg("genomes") = std::vector <Genotype> (), py::arg("low_colour") = 0,
      py::arg("high_colour") = 6);

    m.def("PhenotypeNeighbourhood", &PhenotypeNeighbourhoodAPI,
      "Return the list of the Phenotype_IDs of single mutation neighbours of a genome",
      py::arg("genome"), py::arg("low_colour") = 0, py::arg("high_colour") = 6,
      py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 40,
      py::arg("fixed_table") = true, py::arg("determinism") = 1,
      py::arg("table_file") = std::string("None"));

    m.def("MetricNeighbourhood", &MetricNeighbourhoodAPI,
      "Return metric analysis of the genome neighbourhood as a dictionary",
      py::arg("genome"), py::arg("low_colour") = 0, py::arg("high_colour") = 6,
      py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 40,
      py::arg("fixed_table") = true, py::arg("determinism") = 1,
      py::arg("table_file") = std::string("None"));

    // Metric Sampling

    m.def("MetricSampling", &MetricSamplingAPI,
      "Analysis of the local neighbourhood of mutant genomes isomorphic to the original genome list",
      py::arg("genomes") = std::vector <Genotype> (),
      py::arg("genome_file") = std::string("None"), py::arg("table_file") = std::string("None"),
      py::arg("set_metric_file") = std::string("None"), py::arg("genome_metric_file") = std::string("None"),
      py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 40,
      py::arg("fixed_table") = true, py::arg("determinism") = 1,
      py::arg("n_genes") = 2, py::arg("low_colour") = 0, py::arg("high_colour") = 6,
      py::arg("n_jiggle") = 3, py::arg("dup_aware") = false);

    // Duplicate

    m.def("GenomesDuplication", &GenomesDuplicationAPI,
      "Return a list containing all genomes with a duplicated genes",
      py::arg("genomes") = std::vector <Genotype> ());

    // Phenotype Table

    m.def("PrintTableFromMap", &PrintTableFromMapAPI,
      "Return the GP map built from a list of genomes or a genome file.",
      py::arg("genomes") = std::vector <Genotype> (), py::arg("table_file") = std::string("None"),
      py::arg("gpmap_file") = std::string("None"), py::arg("genome_file") = std::string("None"),
      py::arg("threshold") = 0.25, py::arg("phenotype_builds") = 250,
      py::arg("fixed_table") = true, py::arg("determinism") = 1);

    m.def("LoadTable", &LoadTableAPI,
      "Return a python dictionary of the unordered_map variable of the target PhenotypeTable",
      py::arg("table_file") = std::string("None"));

    // Metric Class

    // py::class_<Genotype_Metrics> (m, "Genotype_Metrics", py::dynamic_attr())
    //   .def(py::init <uint8_t, int8_t, int8_t> ())
    //   .def_readwrite("n_genes", &Genotype_Metrics::n_genes)
    //   .def_readwrite("low_colour", &Genotype_Metrics::low_colour)
    //   .def_readwrite("high_colour", &Genotype_Metrics::high_colour)
    //   .def_readwrite("ref_genotype", &Genotype_Metrics::ref_genotype)
    //   .def_readwrite("original", &Genotype_Metrics::original)
    //   .def_readwrite("ref_pIDs", &Genotype_Metrics::ref_pIDs)
    //   .def_readwrite("pID_counter", &Genotype_Metrics::pID_counter)
    //   .def_readwrite("number_of_neighbours", &Genotype_Metrics::number_of_neighbours)
    //   .def_readwrite("strict_robustness", &Genotype_Metrics::strict_robustness)
    //   .def_readwrite("intersection_robustness", &Genotype_Metrics::intersection_robustness)
    //   .def_readwrite("union_evolvability", &Genotype_Metrics::union_evolvability)
    //   .def_readwrite("complex_diversity", &Genotype_Metrics::complex_diversity)
    //   .def_readwrite("robust_evolvability", &Genotype_Metrics::robust_evolvability)
    //   .def_readwrite("complex_evolvability", &Genotype_Metrics::complex_evolvability)
    //   .def_readwrite("rare", &Genotype_Metrics::rare)
    //   .def_readwrite("unbound", &Genotype_Metrics::unbound)
    //   .def_readwrite("neutral_weight", &Genotype_Metrics::neutral_weight)
    //   .def("to_dict", &Genotype_Metrics::to_dict);
    // m.def("GenerateTable", &GenerateTableAPI, "Test function to print the phenotype Table", py::arg("genotype"), py::arg("filepath") = "None", py::arg("filename") = "None");

    #ifdef VERSION_INFO
      m.attr("__version__") = VERSION_INFO;
    #else
      m.attr("__version__") = "dev";
    #endif
}
