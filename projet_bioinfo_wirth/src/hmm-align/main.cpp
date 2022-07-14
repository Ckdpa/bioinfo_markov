//
// Created by felix on 08/07/22.
//

#include <cstring>
#include "../HMM/HMM.h"

int main(int argc, char *argv[]) {
    HMM hmm(argv[1]);
    Fasta fasta(argv[2]);
    if (argv[3] && !strcmp(argv[3], "--score")) {

    }
    return 0;
}