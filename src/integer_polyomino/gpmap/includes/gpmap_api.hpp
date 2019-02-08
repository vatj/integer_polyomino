#include "integer_model.hpp"
#include "generate.hpp"
#include "metrics.hpp"
#include "io.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// // Assembly Functions
// std::vector <Phenotype_ID> AssemblePlasticGenotypeAPI(Genotype genotype,
//   double threshold, uint16_t phenotype_builds, bool fixed_table,
//   uint8_t determinism, const std::string table_file, py::kwargs kwargs);
//
// std::map <Phenotype_ID, uint16_t> AssemblePlasticGenotypeFrequencyAPI
// (Genotype genotype, double threshold, uint16_t phenotype_builds, bool fixed_table,
//   uint8_t determinism, const std::string table_file, py::kwargs kwargs);
//
// std::vector <std::vector <Phenotype_ID>> AssemblePlasticGenotypesAPI
// (std::vector<Genotype> genomes, double threshold, uint16_t phenotype_builds,
//   bool fixed_table, uint8_t determinism, const std::string table_file, py::kwargs kwargs);
//
// std::vector <std::map <Phenotype_ID, uint16_t>> AssemblePlasticGenotypesFrequencyAPI
// (std::vector<Genotype> genomes, double threshold, uint16_t phenotype_builds,
//   bool fixed_table, uint8_t determinism, const std::string table_file, py::kwargs kwargs);

// Generator Functions
std::vector <Genotype> MinimalGenomesAPI(uint8_t n_genes, int8_t low_colour,
  int8_t high_colour, const std::string genome_file, uint16_t phenotype_builds,
  double threshold, bool fixed_table, uint8_t determinism, py::kwargs kwargs);

void MinimalGenomesVoidAPI(uint8_t n_genes, int8_t low_colour, int8_t high_colour,
  const std::string genome_file, uint16_t phenotype_builds, double threshold,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs);

// GPmap functions
Set_to_Genome MinimalMap(std::vector <Genotype> genomes,
  const std::string table_file, const std::string gpmap_file,
  const std::string genome_file, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs);

String_to_Genome MinimalMapAPI(std::vector <Genotype> genomes,
  const std::string table_file, const std::string gpmap_file,
  const std::string genome_file, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs);

void MinimalMapVoidAPI(std::vector <Genotype> genomes, const std::string table_file,
  const std::string gpmap_file, const std::string genome_file,
  double threshold, uint16_t phenotype_builds, bool fixed_table,
  uint8_t determinism, py::kwargs kwargs);

// Neighbourhood

std::vector<Genotype> GenotypeNeighbourhoodAPI(const Genotype& genome,
  int8_t low_colour, int8_t high_colour, py::kwargs kwargs);

std::vector <std::vector <Phenotype_ID>> PhenotypeNeighbourhoodAPI
(const Genotype genome, int8_t low_colour, int8_t high_colour, double threshold,
  uint16_t phenotype_builds, bool fixed_table, uint8_t determinism,
  const std::string table_file, py::kwargs kwargs);

py::dict MetricNeighbourhoodAPI(Genotype genome,
  int8_t low_colour, int8_t high_colour, double threshold,
  uint16_t phenotype_builds, bool fixed_table, uint8_t determinism,
  const std::string table_file, py::kwargs kwargs);

// Sampling Metrics Functions

void MetricSamplingAPI(std::vector <Genotype> genomes,
  const std::string genome_file, const std::string table_file,
  const std::string set_metric_file, const std::string genome_metric_file,
  double threshold, uint16_t phenotype_builds, bool fixed_table, uint8_t determinism,
  uint8_t n_genes, int8_t low_colour, int8_t high_colour, uint32_t n_jiggle,
  bool dup_aware, py::kwargs kwargs);

// Duplicate Functions

std::vector<Genotype> GenomesDuplicationAPI(std::vector <Genotype> genomes,
py::kwargs kwargs);

// Phenotype Table

void PrintTableFromMapAPI(std::vector <Genotype> genomes,
  const std::string table_file, const std::string gpmap_file,
  const std::string genome_file, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs);

std::unordered_map <uint8_t, std::vector<Phenotype>> LoadTableAPI(const std::string table_file, py::kwargs kwargs);
