//
// Created by felix on 08/07/22.
//


#include "../HMM/HMM.h"

int main(int argc, char *argv[]) {
    // Ne devrait pas arriver, les inputs sont "valides"
    if (argc < 2) {
        return 1;
    }
    HMM hmm(argv[1]);
    hmm.build_print_genseq();
    return 0;
}
