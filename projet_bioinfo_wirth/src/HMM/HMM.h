//
// Created by felix on 09/07/22.
//

#ifndef PROJET_BIOINFO_WIRTH_HMM_H
#define PROJET_BIOINFO_WIRTH_HMM_H


#include <vector>
#include <optional>
#include <map>
#include "Fasta.h"

class HMM {
public:

    /**
     * Ctor from a file and alpha, used for HMM-build
     * @param fasta
     * @param alpha
     */
    explicit HMM(Fasta fasta, float alpha);

    /**
     * Ctor from a file only, used for HMM-genseq
     * @param model_file
     */
    explicit HMM(const std::string& model_file);

    void build_model();

    void build_print_genseq();

    void print_model() const;

private:
    void normalize_matrixes();
    enum class HMMState {
        M = 0,
        D = 1,
        I = 2,
        None = 3
    };
    std::vector<HMMState> build_Pi_k(const std::vector<char>& sequence);
    static std::map<const char, int> alphabet;
    static void display_matrix(std::vector<std::vector<std::optional<float>>> matrix);
    static std::vector<bool> get_marked_columns(const std::vector<std::vector<char>>& sequences, float alpha);
    std::vector<std::vector<char>> sequences_;
    std::vector<bool> marked_columns_;
    std::vector<std::vector<std::optional<float>>> T_;
    std::vector<std::vector<std::optional<float>>> e_M_;
    std::vector<std::vector<std::optional<float>>> e_I_;
    long N_{}; // Nombre d'états
    static char most_probable_char(std::vector<std::optional<float>> &vector);

    static void round_matrix(std::vector<std::vector<std::optional<float>>> matrix);
};


#endif //PROJET_BIOINFO_WIRTH_HMM_H
