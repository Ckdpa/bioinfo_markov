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
    // CTOR
    /**
     * Constructeur depuis un Fasta et un alpha fournis, utilisé pour HMM-build
     * @param fasta le fasta à partir d'une liste de séquences
     * @param alpha paramètre arbitraire entre 0 et 1
     */
    explicit HMM(Fasta fasta, float alpha);

    /**
     * Constructeur à partir d'un fichier .model, utilisé pour HMM-align & HMM-genseq
     * @param model_file
     */
    explicit HMM(const std::string& model_file);

    /**
     * Construction du modèle "Plan 9" selon l'algorithme décrit dans la step 1
     */
    void build_model();

    /**
     * Step 2 - HMM-genseq
     */
    void build_print_genseq();

    /**
     * Écriture du modèle - HMM-build
     */
    void print_model() const;

    /**
     * Algorithme de viterbi & étape retour - HMM-align
     * @param score Vrai s'il faut retourner le score uniquement, faux sinon
     */
    void viterbi(bool score);

    /**
     * Setter
     * @param sequences liste de sequences à déplacer dans la variable sequences_
     */
    void set_sequences(std::vector<std::vector<char>> sequences);

private:
    /**
     * Normalisation des matrices selon les modalités de la step 1
     */
    void normalize_matrixes();
    /**
     * Enum class HMMState
     * Énumère les différents états pouvant apparaître dans la HMM. L'ordre et les valeurs assignées permettent
     * de facilement réaliser des calculs sur les états.
     */
    enum class HMMState {
        M = 0,
        D = 1,
        I = 2,
        None = 3
    };
    // Correspondance globale encodage acide aminé - entier (représenté en size_t pour la cohérence des types)
    static std::map<const char, std::size_t> alphabet;
    // Recherche inversée dans la map ci dessus
    static char find_alphabet_value_of(size_t index);
    // Recherche du caractère d'émission le plus probable du vecteur vector.
    static char most_probable_char(std::vector<std::optional<float>> &vector);
    // Tests sur les résultats non concluants, mais pour la forme
    static void round_matrix(std::vector<std::vector<std::optional<float>>> matrix);
    // HMM-BUILD
    /**
     * Construction de la séquence d'état à partir des séquences de caractères
     * @param sequence liste des séquences à traiter
     * @return un vecteur d'état associé aux séquences
     * @note Utilise la liste de colonnes marquées interne à la classe
     */
    std::vector<HMMState> build_Pi_k(const std::vector<char>& sequence);

    /**
     * Calcule les colonnes marquées et renvoie la liste de booléens associée
     * @param sequences la liste des séquences
     * @param alpha le seuil pour savoir si une colonne est marquée ou non
     * @return une liste de booléen ayant pour valeur true si la colonne est marquée, false sinon
     */
    static std::vector<bool> get_marked_columns(const std::vector<std::vector<char>>& sequences, float alpha);

    // HMM-ALIGN
    /**
     * Écriture de matrice contenant peut-être des flottants
     * @param matrix la matrice à écrire
     */
    static void display_matrix(std::vector<std::vector<std::optional<float>>> matrix);
    /**
     * Surcharge pour une matrice de paire d'enter
     * @param matrix
     */
    [[ maybe_unused ]] static void display_matrix(std::vector<std::vector<std::pair<std::size_t, std::size_t>>> matrix);
    // VARIABLES DE CLASSE
    // Liste de séquences
    std::vector<std::vector<char>> sequences_;
    // Liste de colonnes marquées
    std::vector<bool> marked_columns_;
    // Matrice de probabilités de transition de chaque état vers les suivants à chaque rang de la HMM
    std::vector<std::vector<std::optional<float>>> T_;
    // Matrice de probabilités d'émission de caractère en état M à chaque rang de la HMM
    std::vector<std::vector<std::optional<float>>> e_M_;
    // Matrice de probabilités d'émission de caractère en état I à chaque rang de la HMM
    std::vector<std::vector<std::optional<float>>> e_I_;
    // Nombre d'états de la HMM
    long N_{};
};


#endif //PROJET_BIOINFO_WIRTH_HMM_H
