#include "genotype_phenotype.hpp"

std::vector<Genotype> GenomesDuplication(std::vector<Genotype> genomes);
std::vector<Genotype> GeneDuplication(Genotype& genotype);
std::map<uint8_t,uint8_t> DuplicateGenes(Genotype& genome);
