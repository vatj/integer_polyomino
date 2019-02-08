#include "metrics.hpp"
#include "pybind11/pybind11.h"
#include <sstream>
#include <iterator>
#include <functional>
#include <set>
#include <algorithm>
#include <cstring>

namespace py = pybind11;
using namespace pybind11::literals;

// Main Sampling function

void GP_MapSampler(Set_to_Genome& set_to_genome, PhenotypeTable* pt,
  std::string set_metric_file, std::string genome_metric_file,
  uint8_t n_genes, int8_t low_colour, int8_t high_colour, uint32_t n_jiggle,
  bool dup_aware)
{
  Phenotype_ID unbound_pID = {255, 0}, rare_pID = {0, 0};
  double neutral_weight = 0;
  Genotype mutant;

  uint32_t number_of_genomes = 0;
  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
    number_of_genomes += (iter->second).size();

  py::print("There are", number_of_genomes, "genomes to jiggle a", n_jiggle, "times!");

  // Create new files and open them for writing (erase previous data)
  std::ofstream set_metric_out(set_metric_file);
  std::ofstream genome_metric_out(genome_metric_file);
  header_metric_files(set_metric_out, genome_metric_out);

  for(Set_to_Genome::iterator iter = std::begin(set_to_genome); iter != std::end(set_to_genome); iter++)
  {
    py::print("Currently processing", (iter->second).size(), "genomes for {", "end"_a=" ");
    for(auto pID: iter->first)
      py::print("(", pID.first, ",", pID.second, "),", "end"_a=" ");
    number_of_genomes -= (iter->second).size();
    py::print("}. Only ", number_of_genomes, " left!", "flush"_a=true);

    if(((iter->first).front() == rare_pID) || ((iter->first).back() == unbound_pID))
      continue;

    Set_Metrics set_metrics(n_genes, low_colour, high_colour);
    set_metrics.ref_pIDs = iter->first;

    for(auto genotype: iter->second)
    {
      neutral_weight = ((double) NeutralSize(genotype, 1, high_colour)) / n_jiggle; // This weight will be counted n_jiggle time when added to form the total neutral weight
      set_metrics.originals.emplace_back(genotype);

      #pragma omp parallel for schedule(dynamic) firstprivate(genotype, mutant, neutral_weight)
      for(uint32_t nth_jiggle=0; nth_jiggle<n_jiggle; ++nth_jiggle)
      {
        mutant = JiggleGenotype(genotype, high_colour, dup_aware);

        Genotype_Metrics genome_metric(n_genes, low_colour, high_colour);
        genome_metric.set_reference(mutant, genotype, iter->first, neutral_weight);

        genome_metric.pID_counter = AssemblePlasticGenotypeFrequency(mutant, pt);

        std::vector<Phenotype_ID> pIDs;
        for(auto pID: genome_metric.pID_counter)
          pIDs.emplace_back(pID.first);

        if(pIDs != iter->first)
          set_metrics.misclassified[genotype] = pIDs;

        for(Genotype neighbour : genotype_neighbourhood(mutant, low_colour, high_colour))
        {
           std::vector<Phenotype_ID> neighbour_pIDs = AssemblePlasticGenotype(neighbour, pt);
           genome_metric.analyse_pIDs(neighbour_pIDs);
        }
        #pragma omp critical
        {
           set_metrics.add_genotype_metrics(genome_metric);
        }
      }
    }
    // metrics.emplace_back(set_metrics);
    set_metrics.save_to_file(set_metric_out, genome_metric_out);
  }
  py::print("Metric Sampling has ended!");
}

// Subroutine of the GP_MapSampler

Genotype JiggleGenotype(const Genotype genotype, int8_t high_colour, bool dup_aware)
{
  uint8_t max_colour = high_colour;
  uint8_t min_colour =* std::max_element(genotype.begin(), genotype.end());
  Genotype mutant_genome = genotype;

  if(min_colour + 1 == max_colour)
    return mutant_genome;

  std::vector<uint8_t> neutral_colours( 1 + (max_colour - min_colour) / 2);

  std::generate(neutral_colours.begin() + 1, neutral_colours.end(), [n = min_colour-1] () mutable { return n+=2; });

  std::uniform_int_distribution<size_t> jiggle_index(0, neutral_colours.size() - 1);

  std::map <uint8_t, uint8_t> dups;

  if(dup_aware)
    dups = DuplicateGenes(mutant_genome);

  for(auto& base : mutant_genome)
    base= (base==0) ? neutral_colours[jiggle_index(RNG_Engine)] : base;

  if(dup_aware)
  {
    for(auto dup: dups)
    {
      Genotype::const_iterator first = mutant_genome.begin() + (4 * dup.second);
      Genotype::const_iterator last = mutant_genome.begin() + (4 * (dup.second + 1));
      std::copy(first, last, std::back_inserter(mutant_genome));
    }
  }

  return mutant_genome;
}

std::vector<Genotype> genotype_neighbourhood(const Genotype& genome,
  int8_t low_colour, int8_t high_colour)
{
  std::vector<Genotype> neighbours;
  Genotype neighbour;

  std::vector<int8_t> mutants(high_colour - low_colour + 1);
  std::iota(mutants.begin(), mutants.end(), low_colour);

  for(uint8_t index=0; index<genome.size(); ++index)
  {
    std::swap(*std::find(mutants.begin(), mutants.end(),genome[index]), mutants.back());

    neighbour = genome;
    for(int j=0; j<high_colour - low_colour; ++j)
    {
      neighbour[index] = mutants[j];
      neighbours.emplace_back(neighbour);
    }
  }
  return neighbours;
}

uint64_t nChoosek(uint8_t n, uint8_t k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n - k;
    if (k == 0) return 1;

    uint64_t result = n;

    for(uint8_t i = 2; i <= k; ++i)
    {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

uint64_t combination_with_repetiton(uint8_t space_size , uint8_t sample_size)
{
  if(sample_size==0)
    return 1;
  std::vector<uint8_t> v(sample_size,1);
  uint64_t comb_sum=0;
  while (true)
  {
    for (uint8_t i = 0; i < sample_size; ++i){
      if (v[i] > space_size){
        v[i + 1] += 1;
        for (int16_t k = i; k >= 0; --k)
          v[k] = v[i + 1];
        v[i] = v[i + 1];
      }
    }
    if (v[sample_size] > 0)
      break;
    uint64_t comb_prod=1;
    for(auto x: v)
      comb_prod*=x;
    comb_sum+=comb_prod;
    v[0] += 1;
  }
  return comb_sum;
}

uint64_t NeutralSize(Genotype genotype, uint32_t N_neutral_colours, uint32_t N_possible_interacting_colours)
{
  uint8_t neutral_faces = std::count(genotype.begin(),genotype.end(),0);
  // Clean_Genome(genotype, false);
  std::set<uint8_t> unique_cols(genotype.begin(),genotype.end());

  uint32_t N_interacting_colours= N_possible_interacting_colours, N_interacting_pairs = (unique_cols.size() - 1) / 2;
  uint64_t neutral_interacting=1;

  for(uint8_t n=0; n<N_interacting_pairs; ++n)
    neutral_interacting *= (N_interacting_colours - (2 * n));

  uint32_t N_noninteracting_colours = N_possible_interacting_colours - (unique_cols.size() - 1);
  uint64_t neutral_noninteracting = 0;

  for(uint8_t f=0; f <= neutral_faces; ++f)
  {
    uint64_t pre_sum = nChoosek(neutral_faces, f) * pow(N_neutral_colours, neutral_faces-f);
    uint64_t sum_term=0;

    for(uint8_t U = 0; U <= f; ++U)
    {
      uint64_t pre_prod=1;

      for(uint8_t Un=0;Un<U;++Un)
        pre_prod *= (N_noninteracting_colours - (2 * Un));

      sum_term += pre_prod*combination_with_repetiton(U, f - U);
    }
    neutral_noninteracting += pre_sum * sum_term;
  }
  return neutral_noninteracting*neutral_interacting;
}
