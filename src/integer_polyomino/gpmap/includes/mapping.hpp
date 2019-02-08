#include "core_metrics.hpp"
#include <iostream>
#include <omp.h>
#include "pybind11/pybind11.h"

void PreProcessSampled(std::vector<Genotype> genomes, Set_to_Genome& set_to_genome, PhenotypeTable* pt);
void FilterExhaustive(std::vector<Genotype> genomes, PhenotypeTable* pt);
