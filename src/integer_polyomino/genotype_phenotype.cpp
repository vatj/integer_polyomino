#include "genotype_phenotype.hpp"
#include "pybind11/pybind11.h"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>

namespace py = pybind11;


void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, PhenotypeTable* pt)
{
  Genotype genotype;
  std::vector<Phenotype_ID> pIDs;

  py::print("Mapping", genomes.size(), "genomes, building", pt->phenotype_builds,
    "th times");

  #pragma omp parallel for schedule(dynamic) firstprivate(pIDs, genotype)
  for(uint64_t index=0; index < genomes.size(); index++)
  {
    genotype = genomes[index];
    pIDs = AssemblePlasticGenotype(genotype, pt);

    if(index % 100 ==0)
      py::print("Currently preprocessing genome :", index, "out of", genomes.size());

    #pragma omp critical
      set_to_genome[pIDs].emplace_back(genotype);
  }

  py::print("The GP-Map has been built!");
}

// This function is probably not a priority anymore...
void FilterExhaustive(std::vector<Genotype> genomes, PhenotypeTable* pt)
{
  Genotype genotype;
  std::vector<Genotype> new_genomes;
  std::vector<Phenotype_ID> pIDs;
  Phenotype_ID rare_pID = {0, 0}, unbound_pID = {255, 0};

  std::cout << "Threshold is : " << (ceil(pt->phenotype_builds * pt->UND_threshold));
  std::cout << " out of " <<+ pt->phenotype_builds << " builds \n";

  #pragma omp parallel for schedule(dynamic) firstprivate(pIDs, genotype)
  for(uint64_t index=0; index < genomes.size(); index++)
  {
    genotype = genomes[index];
    pIDs = AssemblePlasticGenotype(genotype, pt);
    if(pIDs.front() == rare_pID || pIDs.back() == unbound_pID)
      continue;
    else
      new_genomes.emplace_back(genotype);

    if(index % 100 ==0)
      std::cout << "Currently filtering genome : " <<+ index << " out of " <<+ genomes.size() << "\n";
  }
  genomes = new_genomes;
}
