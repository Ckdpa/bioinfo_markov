cmake_minimum_required(VERSION 3.18)
project(projet_bioinfo_wirth)

set(CMAKE_CXX_STANDARD 20)

add_executable(HMM-build src/hmm-build/main.cpp src/HMM/Fasta.cpp src/HMM/Fasta.h src/HMM/HMM.cpp src/HMM/HMM.h)
add_executable(HMM-align src/hmm-align/main.cpp src/HMM/Fasta.cpp src/HMM/Fasta.h src/HMM/HMM.cpp src/HMM/HMM.h)
add_executable(HMM-genseq src/hmm-genseq/main.cpp src/HMM/Fasta.cpp src/HMM/Fasta.h src/HMM/HMM.cpp src/HMM/HMM.h)

target_compile_options(HMM-build PRIVATE -g)