#include "genotype_duplicate.hpp"

// Necessary to check that genotype.size() / 4 returns desired value

std::vector<Genotype> GenomesDuplication(std::vector<Genotype> genomes)
{
  std::vector<Genotype> genomes_dup;
  py::print("Adding duplicate genes to", genomes.size(), "genomes");

  for(auto genome: genomes)
    for(auto duplicate: GeneDuplication(genome))
      genomes_dup.emplace_back(duplicate);

  return genomes_dup;
}

std::vector<Genotype> GeneDuplication(Genotype& genotype)
{
  std::vector <Genotype> duplicates;

  for(uint8_t index=0; index < (genotype.size() / 4); ++index)
  {
    Genotype duplicate(genotype.size());
    std::copy(std::begin(genotype), std::end(genotype), std::begin(duplicate));

    for(uint8_t tail=0; tail < 4; tail++)
      duplicate.emplace_back(genotype[(4 * index) + tail]);

    duplicates.emplace_back(duplicate);
  }

  return duplicates;
}

std::map<uint8_t,uint8_t> DuplicateGenes(Genotype& genome) {
  std::map<uint8_t, uint8_t> dups;
  for(int check_index=genome.size()/4-1;check_index>0;--check_index)
    for(int compare_index=0;compare_index<check_index;++compare_index)
      if(std::equal(genome.begin()+check_index*4,genome.begin()+check_index*4+4,genome.begin()+compare_index*4)) {
        genome.erase(genome.begin()+check_index*4,genome.begin()+check_index*4+4);
        dups[check_index]=compare_index;
        break;
      }
  return dups;
}
