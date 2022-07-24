//
// Created by felix on 08/07/22.
//

#include <vector>
#include "../HMM/Fasta.h"
#include "../HMM/HMM.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        return 1;
    }
    HMM hmm(Fasta(argv[1]), std::atof(argv[2]));
    hmm.build_model();
    hmm.print_model();
    return 0;
}


