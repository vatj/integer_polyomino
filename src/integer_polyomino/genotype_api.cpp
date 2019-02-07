#include "genotype_api.hpp"
#include <fstream>

namespace py = pybind11;

// Assembly functions

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

// Minimal Generator functions

std::vector <Genotype> MinimalGenomesAPI(uint8_t n_genes, int8_t low_colour,
  int8_t high_colour, const std::string genome_file, uint16_t phenotype_builds,
  double threshold, bool fixed_table, uint8_t determinism, py::kwargs kwargs)
{
  PhenotypeTable pt;
  PhenotypeTable::UND_threshold = threshold;
  PhenotypeTable::phenotype_builds = phenotype_builds;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";

  std::vector <Genotype> genomes = ExhaustiveMinimalGenotypesFiltered(&pt, n_genes,
     low_colour, high_colour);

  if(genome_file != test)
    PrintGenomeFile(genome_file, genomes);

  return genomes;
}

void MinimalGenomesVoidAPI(uint8_t n_genes, int8_t low_colour, int8_t high_colour,
  const std::string genome_file, uint16_t phenotype_builds, double threshold,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs)
{
  std::vector <Genotype> genomes = MinimalGenomesAPI(n_genes, low_colour,
    high_colour, genome_file, phenotype_builds, threshold, fixed_table,
    determinism, kwargs);
}


// Minimal GP_map

Set_to_Genome MinimalMap(std::vector <Genotype> genomes,
  const std::string table_file, const std::string gpmap_file,
  const std::string genome_file, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs)
{
  PhenotypeTable pt;
  PhenotypeTable::UND_threshold = threshold;
  PhenotypeTable::phenotype_builds = phenotype_builds;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";

  Set_to_Genome gpmap;

  if(genome_file != test)
    LoadGenomeFile(genome_file, genomes);

  if(table_file != test)
    pt.LoadTable(table_file);

  PreProcessSampled(genomes, gpmap, &pt);

  if(gpmap_file != test)
    PrintPreProcessFile2(gpmap_file, gpmap);

  return gpmap;
}

String_to_Genome MinimalMapAPI(std::vector <Genotype> genomes,
  const std::string table_file, const std::string gpmap_file,
  const std::string genome_file, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs)
{
  Set_to_Genome gpmap = MinimalMap(genomes, table_file, gpmap_file,
    genome_file, threshold, phenotype_builds, fixed_table, determinism, kwargs);

  String_to_Genome gpmap2;
  std::stringstream pIDs;

  for(auto kv: gpmap)
  {

    pIDs << "{";
    for(auto pID: kv.first)
      pIDs << "(" <<+ pID.first << "," <<+ pID.second << "),";

    gpmap2[pIDs.str().substr(0, pIDs.str().size() - 1) + "}"] = kv.second;
    pIDs.str(std::string());
  }

  return gpmap2;
}

void MinimalMapVoidAPI(std::vector <Genotype> genomes, const std::string table_file,
  const std::string gpmap_file, const std::string genome_file, double threshold,
  uint16_t phenotype_builds, bool fixed_table, uint8_t determinism, py::kwargs kwargs)
{
  Set_to_Genome gpmap = MinimalMap(genomes, table_file, gpmap_file,
    genome_file, threshold, phenotype_builds, fixed_table, determinism, kwargs);
}

// Local neighbourhood

std::vector<Genotype> GenotypeNeighbourhoodAPI(const Genotype& genome,
  int8_t low_colour, int8_t high_colour, py::kwargs kwargs)
{
  return genotype_neighbourhood(genome, low_colour, high_colour);
}

std::vector <std::vector <Phenotype_ID>> PhenotypeNeighbourhoodAPI
(const Genotype genome, int8_t low_colour, int8_t high_colour, double threshold,
  uint16_t phenotype_builds, bool fixed_table, uint8_t determinism,
  const std::string table_file, py::kwargs kwargs)
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

  for(auto genome: genotype_neighbourhood(genome, low_colour, high_colour))
    vector_pIDs.emplace_back(AssemblePlasticGenotype(genome, &pt));

  return vector_pIDs;
}

py::dict MetricNeighbourhoodAPI(Genotype genome,
  int8_t low_colour, int8_t high_colour, double threshold,
  uint16_t phenotype_builds, bool fixed_table, uint8_t determinism,
  const std::string table_file, py::kwargs kwargs)
{
  PhenotypeTable pt = PhenotypeTable();
  PhenotypeTable::phenotype_builds = phenotype_builds;
  PhenotypeTable::UND_threshold = threshold;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";
  std::vector <Phenotype_ID> pIDs;
  Genotype_Metrics genome_metric(genome.size() / 4, low_colour, high_colour);

  if(table_file != test)
  {
    py::print("Loading table from : ", table_file);
    pt.LoadTable(table_file);
  }

  genome_metric.pID_counter = AssemblePlasticGenotypeFrequency(genome, &pt);
  std::vector <Phenotype_ID> ref_pIDs;

  for(auto pID: genome_metric.pID_counter)
    ref_pIDs.emplace_back(pID.first);

  genome_metric.set_reference(genome, genome, ref_pIDs, NeutralSize(genome, 1, high_colour));

  for(auto genotype: genotype_neighbourhood(genome, low_colour, high_colour)) {
    pIDs = AssemblePlasticGenotype(genotype, &pt);
    genome_metric.analyse_pIDs(pIDs);
  }

  return genome_metric.to_dict();
}


// Metric Sampling

void MetricSamplingAPI(std::vector <Genotype> genomes,
  const std::string genome_file, const std::string table_file,
  const std::string set_metric_file, const std::string genome_metric_file,
  double threshold, uint16_t phenotype_builds, bool fixed_table, uint8_t determinism,
  uint8_t n_genes, int8_t low_colour, int8_t high_colour, uint32_t n_jiggle,
  bool dup_aware, py::kwargs kwargs)
{
  PhenotypeTable pt;
  PhenotypeTable::UND_threshold = threshold;
  PhenotypeTable::phenotype_builds = phenotype_builds;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";
  Set_to_Genome gpmap;

  if(genome_file != test)
    LoadGenomeFile(genome_file, genomes);

  if(table_file != test)
    pt.LoadTable(table_file);

  PreProcessSampled(genomes, gpmap, &pt);

  GP_MapSampler(gpmap, &pt, set_metric_file, genome_metric_file,
    n_genes, low_colour, high_colour, n_jiggle, dup_aware);
}

// Duplicate genes

std::vector<Genotype> GenomesDuplicationAPI(std::vector <Genotype> genomes,
  py::kwargs kwargs)
{
  return GenomesDuplication(genomes);
}

// Phenotype Table

void PrintTableFromMapAPI(std::vector <Genotype> genomes,
  const std::string table_file, const std::string gpmap_file,
  const std::string genome_file, double threshold, uint16_t phenotype_builds,
  bool fixed_table, uint8_t determinism, py::kwargs kwargs)
{
  PhenotypeTable pt;
  PhenotypeTable::UND_threshold = threshold;
  PhenotypeTable::phenotype_builds = phenotype_builds;
  pt.FIXED_TABLE = fixed_table;
  Phenotype::DETERMINISM_LEVEL = determinism;
  std::string test = "None";

  Set_to_Genome gpmap;

  if(genome_file != test)
    LoadGenomeFile(genome_file, genomes);

  PreProcessSampled(genomes, gpmap, &pt);

  if(table_file != test) {
    py::print("Printing table to :", table_file);
    pt.PrintTable(table_file);
  }

  if(gpmap_file != test)
    PrintPreProcessFile2(gpmap_file, gpmap);
}

std::unordered_map <uint8_t, std::vector<Phenotype>> LoadTableAPI(const std::string table_file,
  py::kwargs kwargs)
{
  PhenotypeTable pt;
  pt.LoadTable(table_file);

  return pt.known_phenotypes;
}
