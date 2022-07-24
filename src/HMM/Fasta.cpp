//
// Created by felix on 09/07/22.
//

#include "Fasta.h"

Fasta::Fasta(const std::string& filename):
file_(filename){
}

std::vector<std::vector<char>> Fasta::parse() {
    std::vector<std::vector<char>> sequences;
    std::string line;
    std::vector<char> sequence;
    while (getline(this->file_, line)) {
        if (line.find('>') != std::string::npos) {
            if (not sequence.empty()) {
                sequences.emplace_back(sequence);
            }
            sequence.clear();
        }
        else {
            copy(line.begin(), line.end(), back_inserter(sequence));
        }
    }
    if (not sequence.empty()) {
        sequences.emplace_back(sequence);
    }
    return sequences;
}
