//
// Created by felix on 09/07/22.
//

#ifndef PROJET_BIOINFO_WIRTH_FASTA_H
#define PROJET_BIOINFO_WIRTH_FASTA_H


#include <string>
#include <fstream>
#include <vector>

class Fasta {
public:
    explicit Fasta(const std::string& filename);
    std::vector<std::vector<char>> parse();
private:
    std::ifstream file_;
};


#endif //PROJET_BIOINFO_WIRTH_FASTA_H
