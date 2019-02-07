// #include "genotype_generate.hpp"
#include "genotype_core_metrics.hpp"
#include "genotype_iofunc.hpp"
#include <iostream>

// namespace simulation_params
// {
//   extern uint16_t n_genes, colours, metric_colours;
//   extern uint32_t n_jiggle;
//   extern std::mt19937 RNG_Engine;
//   extern bool dup_aware;
// }

// namespace io_params
// {
//   extern std::string set_metric_file, genome_metric_file;
// }

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
