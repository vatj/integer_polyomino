#include "core_metrics.hpp"
#include "io.hpp"
#include <iostream>

void GP_MapSampler(Set_to_Genome& set_to_genome, PhenotypeTable* pt,
  std::string set_metric_file, std::string genome_metric_file,
  uint8_t n_genes, int8_t low_colour, int8_t high_colour, uint32_t n_jiggle,
  bool dup_aware);

// Subroutines of the metric sampler
std::vector<Genotype> genotype_neighbourhood(const Genotype& genome,
  int8_t low_colour, int8_t high_colour);
Genotype JiggleGenotype(const Genotype genotype, int8_t high_colour, bool dup_aware);

/*Neutral size calculations*/
uint64_t NeutralSize(Genotype genotype,uint32_t N_neutral_colours,uint32_t N_possible_interacting_colours);
uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size);
uint64_t nChoosek(uint8_t n, uint8_t k);
