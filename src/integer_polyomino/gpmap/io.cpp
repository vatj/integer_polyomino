#include "io.hpp"
#include "pybind11/pybind11.h"
#include <sstream>
#include <iterator>
#include <string>
#include <utility>

namespace py = pybind11;


void PrintGenomeFile(std::string genome_file, std::vector<Genotype>& genomes)
{
  py::print("Printing", genomes.size(), "genomes to file :", genome_file);
  std::ofstream fout(genome_file);
  for(auto genome: genomes)
  {
    for(auto base: genome)
      fout <<+ base << " ";
    fout << "\n";
  }
}

void PrintPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome)
{
  py::print("Printing preprocessed genomes to file :", preprocess_file);
  std::ofstream fout(preprocess_file);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    fout << "x ";
    for (auto pID: iter->first)
      fout <<+ pID.first << " " <<+ pID.second << " ";
    fout << "\n";

    for(auto genome: iter->second)
    {
      for(auto index: genome)
        fout <<+ index << " ";
      fout << "\n";
    }
  }
}

void PrintPreProcessFile2(std::string preprocess_file, Set_to_Genome& set_to_genome)
{
  py::print("Printing preprocessed genomes to file :", preprocess_file);
  std::ofstream fout(preprocess_file);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    for(auto genome: iter->second)
    {
      for(auto index: genome)
        fout <<+ index << " ";
      fout << "{";
      for (auto pID: iter->first)
        fout <<+ pID.first << " " <<+ pID.second << " ";
      fout << "}\n";
    }
  }
}

void PrintSetTable(std::string set_file, Set_to_Genome& set_to_genome)
{
  py::print("Printing set table to file :", set_file);
  std::ofstream fout(set_file);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    fout << "{";
    for (auto pID: iter->first)
      fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
    fout.seekp((long) fout.tellp() - 1);
    fout << "} ";

    fout << "[";
    for (auto original: iter->second)
    {
      fout << "(";
      for (auto face: original)
        fout <<+ face << ",";
      fout.seekp((long) fout.tellp() - 1);
      fout << "),";
    }
    fout.seekp((long) fout.tellp() - 1);
    fout << "]\n";
  }
}

void header_metric_files(std::ofstream& set_metric_out, std::ofstream& genome_metric_out)
{
  // Logging
  // py::print("Print metrics to files :");
  // py::print(set_metric_file);
   // py::print(genome_metric_file);

  // Header for the metric files
  set_metric_out << "srobustness irobustness evolvability";
  set_metric_out << " robust_evolvability complex_evolvability";
  set_metric_out << " complex_diversity rare unbound";
  set_metric_out << " analysed misclassified neutral_size";
  set_metric_out << " diversity diversity_tracker originals misclassified_details pIDs\n";

  genome_metric_out << "genome original srobustness irobustness";
  genome_metric_out << " evolvability complex_evolvability robust_evolvability";
  genome_metric_out << " rare unbound complex_diversity diversity";
  genome_metric_out << " neutral_weight frequencies pIDs\n";
}

void LoadGenomeFile(std::string genome_file, std::vector<Genotype>& genomes)
{
  std::string str;
  Genotype genotype;
  std::ifstream genome_in(genome_file);

  py::print("Loading genome from file :", genome_file);

  while (std::getline(genome_in, str))
  {
    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    genomes.emplace_back(genotype);
  }
}


void LoadPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome)
{
  std::string str;
  Genotype genotype;
  std::vector<int> pre_pIDs;
  std::vector<Phenotype_ID> pIDs;
  std::ifstream fin(preprocess_file);

  py::print("Loading preprocess file :", preprocess_file);

  while (std::getline(fin, str))
  {
    if(str.compare(0, 1, "x"))
    {
      std::istringstream is(str);
      pre_pIDs.assign(std::istream_iterator<int>(is), std::istream_iterator<int>());

      for(uint8_t index=0; index < pre_pIDs.size() - 1; index+=2)
        pIDs.emplace_back(std::make_pair(pre_pIDs[index], pre_pIDs[index + 1]));
      continue;
    }

    std::istringstream is( str );
    genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    set_to_genome[pIDs].emplace_back(genotype);
  }
}

void LoadPhenotypeTable(std::string phenotype_file, PhenotypeTable* pt_it)
{
  std::ifstream pheno_in(phenotype_file);
  std::string str;

  py::print("Loading phenotype table :", phenotype_file);

  while (std::getline(pheno_in, str))
  {
    std::stringstream iss(str);
    int number;
    std::vector<uint8_t> phenotype_line;
    while (iss>>number)
      phenotype_line.push_back(static_cast<uint8_t>(number));
    Phenotype phen;
    phen.dx=phenotype_line[2];
    phen.dy=phenotype_line[3];
    phen.tiling=std::vector<uint8_t>(phenotype_line.begin()+4, phenotype_line.end());
    pt_it->known_phenotypes[phenotype_line[0]].emplace_back(phen);
  }
}
