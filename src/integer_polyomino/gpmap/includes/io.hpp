#include "integer_model.hpp"
#include "duplicate.hpp"
#include <iostream>


void PrintGenomeFile(std::string genome_file, std::vector<Genotype>& genomes);
void PrintPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome);
void PrintSetTable(std::string set_file, Set_to_Genome& set_to_genome);
void PrintPreProcessFile2(std::string preprocess_file, Set_to_Genome& set_to_genome);

void header_metric_files(std::ofstream& set_metric_out, std::ofstream& genome_metric_out);

void LoadGenomeFile(std::string genome_file, std::vector<Genotype>& genomes);
void LoadPreProcessFile(std::string preprocess_file, Set_to_Genome& set_to_genome);
void LoadPhenotypeTable(std::string phenotype_file, PhenotypeTable* pt_it);
