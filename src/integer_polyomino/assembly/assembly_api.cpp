#include "assembly_api.hpp"

namespace py = pybind11;


std::vector <Phenotype_ID> AssemblePlasticGenotypeAPI(Genotype genotype,
  double threshold, uint16_t phenotype_builds, bool fixed_table,
  uint8_t determinism, const std::string table_file, py::kwargs kwargs)
{
  PhenotypeTable pt = PhenotypeTable();
  PhenotypeTable::phenotype_builds = phenotype_builds;
  PhenotypeTable::UND_threshold = threshold;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";

  if(table_file != test)
  {
    py::print("Loading table from : ", table_file);
    pt.LoadTable(table_file);
  }

  return AssemblePlasticGenotype(genotype, &pt);
}

std::map <Phenotype_ID, uint16_t> AssemblePlasticGenotypeFrequencyAPI(Genotype genotype,
  double threshold, uint16_t phenotype_builds, bool fixed_table,
  uint8_t determinism, const std::string table_file, py::kwargs kwargs)
{
  PhenotypeTable pt = PhenotypeTable();
  PhenotypeTable::phenotype_builds = phenotype_builds;
  PhenotypeTable::UND_threshold = threshold;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";

  if(table_file != test)
  {
    py::print("Loading table from : ", table_file);
    pt.LoadTable(table_file);
  }

  return AssemblePlasticGenotypeFrequency(genotype, &pt);
}

std::vector <std::vector <Phenotype_ID>> AssemblePlasticGenotypesAPI
(std::vector<Genotype> genomes, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, const std::string table_file, py::kwargs kwargs)
{
  PhenotypeTable pt = PhenotypeTable();
  PhenotypeTable::phenotype_builds = phenotype_builds;
  PhenotypeTable::UND_threshold = threshold;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";
  std::vector <std::vector <Phenotype_ID>> vector_pIDs;

  if(table_file != test)
  {
    py::print("Loading table from : ", table_file);
    pt.LoadTable(table_file);
  }

  for(auto genome: genomes)
    vector_pIDs.emplace_back(AssemblePlasticGenotype(genome, &pt));

  return vector_pIDs;
}

std::vector <std::map <Phenotype_ID, uint16_t>> AssemblePlasticGenotypesFrequencyAPI
(std::vector<Genotype> genomes, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, const std::string table_file, py::kwargs kwargs)
{
  PhenotypeTable pt = PhenotypeTable();
  PhenotypeTable::phenotype_builds = phenotype_builds;
  PhenotypeTable::UND_threshold = threshold;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";
  std::vector <std::map <Phenotype_ID, uint16_t>> vector_pIDcounter;

  if(table_file != test)
  {
    py::print("Loading table from : ", table_file);
    pt.LoadTable(table_file);
  }

  for(auto genome: genomes)
    vector_pIDcounter.emplace_back(AssemblePlasticGenotypeFrequency(genome, &pt));

  return vector_pIDcounter;
}
