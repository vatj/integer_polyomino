#include "integer_model.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// Assembly Functions
std::vector <Phenotype_ID> AssemblePlasticGenotypeAPI(Genotype genotype,
  double threshold, uint16_t phenotype_builds, bool fixed_table,
  uint8_t determinism, const std::string table_file, py::kwargs kwargs);

std::map <Phenotype_ID, uint16_t> AssemblePlasticGenotypeFrequencyAPI
(Genotype genotype, double threshold, uint16_t phenotype_builds, bool fixed_table,
  uint8_t determinism, const std::string table_file, py::kwargs kwargs);

std::vector <std::vector <Phenotype_ID>> AssemblePlasticGenotypesAPI
(std::vector<Genotype> genomes, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, const std::string table_file, py::kwargs kwargs);

std::vector <std::map <Phenotype_ID, uint16_t>> AssemblePlasticGenotypesFrequencyAPI
(std::vector<Genotype> genomes, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, const std::string table_file, py::kwargs kwargs);
