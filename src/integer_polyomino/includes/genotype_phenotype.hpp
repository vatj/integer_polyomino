#include "genotype_core_metrics.hpp"
#include <iostream>
#include <omp.h>
#include "pybind11/pybind11.h"

std::vector<Phenotype_ID> GetSetPIDs(Genotype genotype, PhenotypeTable* pt_it);
std::map<Phenotype_ID, uint16_t> GetPIDCounter(Genotype genotype, PhenotypeTable* pt_it);

void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, PhenotypeTable* pt);
void FilterExhaustive(std::vector<Genotype> genomes, PhenotypeTable* pt);
