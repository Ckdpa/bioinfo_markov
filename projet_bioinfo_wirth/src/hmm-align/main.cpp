//
// Created by felix on 08/07/22.
//

#include <cstring>
#include <iostream>
#include "../HMM/HMM.h"

int main(int argc, char *argv[]) {
    bool score = false;
    int begin = 1;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--score")) {
            score = true;
            begin = i + 1;
        }
    }
    HMM hmm(argv[begin]);
    hmm.set_sequences(Fasta(argv[begin + 1]).parse());
    hmm.viterbi(score);
    return 0;
}