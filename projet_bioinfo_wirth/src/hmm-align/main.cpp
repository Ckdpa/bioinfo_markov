//
// Created by felix on 08/07/22.
//

#include <cstring>
#include <iostream>
#include "../HMM/HMM.h"

int main(int argc, char *argv[]) {
    HMM hmm(argv[1]);
    hmm.set_sequences(Fasta(argv[2]).parse());
    std::cout << "VITERBI\n";
    hmm.viterbi();
    std::cout << "VITERBOOK\n";
    hmm.viterbook();
    if (argv[3] && !strcmp(argv[3], "--score")) {

    }
    return 0;
}