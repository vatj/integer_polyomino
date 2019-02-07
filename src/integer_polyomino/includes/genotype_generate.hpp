#include "genotype_iofunc.hpp"
#include <iostream>

// namespace simulation_params
// {
//   extern uint16_t n_genes, colours;
//   extern std::mt19937 RNG_Engine;
//   extern uint16_t phenotype_builds;
//   extern uint32_t n_samples;
//   extern bool allow_duplicates;
// }


std::vector<Genotype> ExhaustiveMinimalGenotypesFiltered(PhenotypeTable* pt,
  uint8_t n_genes, int8_t low_colour, int8_t high_colour);

std::vector<Genotype> ExhaustiveMinimalGenotypesIL(PhenotypeTable* pt);
std::vector<Genotype> SampleMinimalGenotypes(PhenotypeTable* pt);
std::vector<Genotype> ExhaustiveMinimalGenotypesFastFiltered(PhenotypeTable* pt);
// std::vector<Genotype> ExhaustiveFullGenotypes2(uint16_t colours, PhenotypeTable* pt);

// Temporary Fix
std::vector<Genotype> ExhaustiveMinimalGenotypesFilteredDuplicate(std::vector<Genotype>& genomes, PhenotypeTable* pt);


/*Minimal genotype methods*/
struct NecklaceFactory
{
  int8_t low_colours,high_colours=1;
  std::vector<std::vector<int8_t> > necklaces;
  std::vector<int8_t> necklace_grower;

  NecklaceFactory(int8_t lcol,int8_t hcol) {
    low_colours=lcol;
    high_colours=hcol;
    necklace_grower.assign(5,low_colours);
    crsms_gen(1,1);
  }

  bool is_finite_necklace(std::vector<int8_t>& neck);

  void is_necklace(int64_t j);

  void crsms_gen(int64_t n, int64_t j);
};

struct GenotypeGenerator
{
  bool is_done=false;
  uint8_t n_genes;
  int8_t low_colours,high_colours;
  std::vector<int32_t> necklace_states;
  std::vector<std::vector<int8_t> > necklaces;
  uint32_t n_necklaces;


  GenotypeGenerator(uint8_t a, int8_t b, int8_t c) {
    n_genes= a;
    low_colours=b;
    high_colours=c;
    necklace_states.assign(a,0);
    NecklaceFactory necks=NecklaceFactory(low_colours,high_colours);
    necklaces=necks.necklaces;
    n_necklaces=necklaces.size();
    necklace_states[0]=1;
  }

  Genotype operator() ()
  {
    return !is_done ? next_genotype() : Genotype{};
  }

  bool valid_growing_faces(Genotype& genotype, int8_t max_face)
  {
    for(auto face : genotype)
    {
      if(face>max_face)
        return false;
      if(face==max_face)
        max_face+=2;
    }
    return true;
  }

  bool valid_bindings(Genotype& genotype)
  {
    //negative interfaces self-interact, trivially valid bindings
    for(int8_t interface=1; interface<=*std::max_element(genotype.begin(), genotype.end()); interface+=2)
    {
      if(std::find(genotype.begin(),genotype.end(),interface)!=genotype.end())
      { //is present
        if(std::find(genotype.begin(),genotype.end(),interface+1)==genotype.end()) //is not present
          return false;
      }
      else {
        if(std::find(genotype.begin(),genotype.end(),interface+1)!=genotype.end())
          return false;
      }
    }
    return true;
  }

  bool valid_genotype(Genotype& genotype)
  {
    return (valid_growing_faces(genotype,1) && valid_bindings(genotype));
  }

  void increment_states(std::vector<int32_t>& states)
  {
    ++states.back();
    int32_t zero_state_init=states[0];
    for(int32_t rind=states.size()-1;rind>=0;--rind) {
      if(states[rind]>=static_cast<int32_t>(n_necklaces)) {
        if(rind==0) {
          is_done=true;
          return;
        }
        else {
          states[rind]=0;
          ++states[rind-1];
        }
      }
    }
    /*! short cuts based on 0001, 0013 etc ideas, not valid in negative number
      would need new offsets to determine location in necklace chain*/
    if(low_colours==0) {
      if(zero_state_init==1 && states[0]!=1)
        states[1]=high_colours+3;
      if(zero_state_init==(high_colours+3) && states[0]!=(high_colours+3)) {
        is_done=true;
        return;
      }
    }
    auto max_iter=std::max_element(states.begin(),states.end());
    std::replace(max_iter,states.end(),0,*max_iter);
  }

  Genotype next_genotype() {
    Genotype genotype;

    while(!is_done) {
    inc_lab:
      genotype.clear();
      genotype.reserve(n_genes*4);

      for(uint8_t index=0; index!= n_genes;++index) {
        int8_t input_face=*std::max_element(genotype.begin(),genotype.end());
        int8_t next_max_face= input_face>0 ? ((input_face-1)/2)*2+3 : 1;

        while(!valid_growing_faces(necklaces[necklace_states[index]], next_max_face)) {
          const uint32_t post_inc = necklace_states[index]+1;
          if(post_inc==necklaces.size()) {
            increment_states(necklace_states);
            goto inc_lab;
          }
          std::fill(necklace_states.begin()+index,necklace_states.end(),post_inc);
        }

        genotype.insert(genotype.end(),necklaces[necklace_states[index]].begin(),necklaces[necklace_states[index]].end());
      }

      increment_states(necklace_states);

      if(valid_genotype(genotype))
        return genotype;
    }
    genotype.clear();
    return genotype;
  }
};
