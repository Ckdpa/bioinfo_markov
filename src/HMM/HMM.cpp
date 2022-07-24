//
// Created by felix on 09/07/22.
//

#include <iostream>
#include <iomanip>
#include <utility>
#include <valarray>
#include "HMM.h"

std::map<const char, std::size_t> HMM::alphabet = {{'A', 0}, {'C', 1}, {'D', 2}, {'E', 3},
                                                        {'F', 4}, {'G', 5}, {'H', 6}, {'I', 7},
                                                        {'K', 8}, {'L', 9}, {'M', 10}, {'N', 11},
                                                        {'P', 12}, {'Q', 13}, {'R', 14},
                                                        {'S', 15}, {'T', 16}, {'V', 17},
                                                        {'W', 18}, {'Y', 19}};


// UTILS
// Arrondir utilisé (avant d'écrire les matrices)
inline float round(float val )
{
    if( val < 0 ) return ceil(val - 0.5);
    return floor(val + 0.5);
}

// Trouver le maximum d'un vecteur entre un point de départ et d'arrivée (ou tout le vecteur par défaut.)
// Le i_factor permet d'appliquer une puissance toutes les valeurs = 2 [mod 3]. Utile pour HMM-genseq.
// Les paramètres par défaut permettent de l'utiliser sans prêter attention à ces particularités, auquel cas on retourne
// simplement le maximum du vecteur
std::size_t index_of_max(std::vector<std::optional<float>> & vector, int i_factor=1, std::size_t start=0, std::size_t stop=0);
// Implémentation
std::size_t index_of_max(std::vector<std::optional<float>> & vector, int i_factor, std::size_t start, std::size_t stop) {

    auto value = 0.;
    if (start == stop) {
        start = 0;
        stop = vector.size();
    }
    std::size_t max_index = start;
    for (auto index = start; index < stop; index++) {
        // Raise to the power i_factor (used when searching for the next state)
        if (index % 3 == 2 && i_factor != 1) {
            if (std::pow(vector[index].value(), i_factor) > value) {
                value = std::pow(vector[index].value(), i_factor);
                max_index = index;
            }
        }
        else if (vector[index].value() > value) {
            value = vector[index].value();
            max_index = index;
        }
    }
    return max_index;
}

// Constructeur à partir d'un fichier fasta et d'un alpha (utilisé pour HMM-build)
HMM::HMM(Fasta fasta, float alpha)
:sequences_(fasta.parse()),
marked_columns_(get_marked_columns(sequences_, alpha)),
N_(std::count(marked_columns_.begin(),marked_columns_.end(), true) + 1 )
{
    T_.emplace_back(9, 1.);
    // Le premier état n'émet pas de caractère, aussi la ligne est remplie de NaN (ici, optional sans valeur)
    e_M_.emplace_back(20, std::optional<float>());
    e_I_.emplace_back(20, 1.);
    for (auto i = 0; i < N_-2; i++) {
        T_.emplace_back(9, 1.);
        e_M_.emplace_back(20, 1.);
        e_I_.emplace_back(20, 1.);
    }
    // Dernière ligne : 0 pour T
    T_.emplace_back(9, 0.);
    e_M_.emplace_back(20, 1.);
    e_I_.emplace_back(20, 1.);
    // L'état D0 n'existe pas, aussi les transitions à partir de cet état sont mises à 0.
    T_[0][3] = 0.;
    T_[0][4] = 0.;
    T_[0][5] = 0.;
}

// Constructeur à partir d'un fichier (utilisé par HMM-genseq, HMM-align)
HMM::HMM(const std::string& model_file) {
    // parse the model
    std::ifstream input(model_file);
    input >> N_;
    // Initialisation des matrices à 0, première colonne séparément pour gérer e_M[0]
    T_.emplace_back(9, 0.);
    e_M_.emplace_back(20, std::optional<float>());
    e_I_.emplace_back(20, 0.);
    for (auto i = 0; i < N_-1; i++) {
        T_.emplace_back(9, 0.);
        e_M_.emplace_back(20, 0.);
        e_I_.emplace_back(20, 0.);
    }
    std::string test;
    //T
    // read line by line
    for (auto model_line = 0; model_line < N_; model_line++) {
        // Get the next line
        input >> test;
        for (auto model_column = 0; model_column < 9; model_column++) {
            T_[model_line][model_column] = std::stof(test.substr(6 * model_column, 5));
        }
    }
    //skip NAN line
    input >> test;
    //e_M
    // read line by line
    for (auto model_line = 1; model_line < N_; model_line++) {
        // Get the next line
        input >> test;
        for (auto model_column = 0; model_column < 20; model_column++) {
            e_M_[model_line][model_column] = std::stof(test.substr(6 * model_column, 5));
        }
    }
    //e_i
    // read line by line
    for (auto model_line = 0; model_line < N_; model_line++) {
        // Get the next line
        input >> test;
        for (auto model_column = 0; model_column < 20; model_column++) {
            e_I_[model_line][model_column] = std::stof(test.substr(6 * model_column, 5));
        }
    }
}

// Calcul des colonnes marquées : Si une colonne a en moyenne moins de alpha '-', elle est marquée.
std::vector<bool> HMM::get_marked_columns(const std::vector<std::vector<char>>& sequences, float alpha) {
    std::vector<bool> ret;
    std::vector<float> acc(sequences[0].size(), 0.0);
    for (auto seq : sequences) {
        for (std::size_t i = 0; i < seq.size(); i++){
            if (seq[i] == '-') {
                acc[i]++;
            }
        }
    }
    for (auto count : acc) {
        if (count / static_cast<float>(sequences.size()) < alpha) {
            ret.push_back(true);
        } else {
            ret.push_back(false);
        }
    }
    return ret;
}

// Écriture du modèle comme demandé dans HMM-build
void HMM::print_model() const {
#ifdef DEBUG
    for (bool mark : marked_columns_) {
        std::cout << mark << " ";
    }
    std::cout<< std::endl;
#endif
    // 3 chiffres après la virgule
    std::cout << std::setprecision(3) << std::fixed;
    std::cout << N_ << std::endl;
    display_matrix(T_);
    display_matrix(e_M_);
    display_matrix(e_I_);
}

// Écriture de matrice, ligne apr ligne, si la valeur existe
void HMM::display_matrix(std::vector<std::vector<std::optional<float>>> matrix) {
    for (auto & line : matrix) {
        for (std::size_t j = 0; j < line.size(); j++) {
            if (line[j].has_value()) {
                std::cout << line[j].value();
            } else {
                std::cout << "nan";
            }
            if (j != line.size() - 1) {
                std::cout << ',';
            } else {
                std::cout << std::endl;
            }
        }
    }
}

// Algorithme de construction de modèle (HMM-build)
void HMM::build_model() {
    // Pour chaque séquence Ak de A
    std::vector<HMMState> Pi_k;
    std::size_t l_count; // l
    int model_column; // u
    // Pour chaque séquence de A_k
    for (auto & A_k : sequences_) {
        // Début au rang 0
        model_column = 0;
        // Détermination de Pi_k
        Pi_k.clear();
        Pi_k = build_Pi_k(A_k);
        // Fin détermination de Pi_k

        // Mettre à jour les matrices à partir de Pi_k
        l_count = 0;
#ifdef DEBUG
        for (auto& state : Pi_k) {
            std::cout << static_cast<int>(state) << " ";
        }
        std::cout << std::endl;
#endif
        // Déterminer la première position l0 non None
        while (l_count < Pi_k.size() && Pi_k[l_count] == HMMState::None) {
            l_count++;
        }

        // Si le premier état l0 non None est I
        if (Pi_k[l_count] == HMMState::I) {
            // Ajouter 1 à la position correspondant à l'acide aminé A_k[l0] dans e_I[0]
            e_I_[model_column][alphabet[A_k[l_count]]] =
                    e_I_[model_column][alphabet[A_k[l_count]]].value() + 1;
        }
        // Mettre à jour T0 en considérant que l'état précédent est M
        T_[model_column][static_cast<int>(Pi_k[l_count])] =
                T_[model_column][static_cast<int>(Pi_k[l_count])].value() + 1;
        // Ensuite, pour chaque colonne l >= l0 de A_k
        while (l_count < A_k.size()) {
            if (marked_columns_[l_count]) {
                model_column++;
            }
            // Si Pi_k[l] est un état M
            // Et Ak_l différent de '-'
            if (Pi_k[l_count] == HMMState::M && A_k[l_count] != '-' && model_column != 0) {
                // Ajouter 1 à la position correspondant à l'acide aminé A_k[l] dans e_M[u]
                // Skip caractère "X"
                if (A_k[l_count] != 'X') {
                    e_M_[model_column][alphabet[A_k[l_count]]] =
                            e_M_[model_column][alphabet[A_k[l_count]]].value() + 1;
                }
            }
            // Respectivement I
            else if (Pi_k[l_count] == HMMState::I && A_k[l_count] != '-') {
                // Et Ak_l différent de '-'
                // Ajouter 1 à la position correspondant à l'acide aminé A_k[l] dans e_I[u]
                // Skip le caractère "X"
                if (A_k[l_count] != 'X') {
                    e_I_[model_column][alphabet[A_k[l_count]]] =
                            e_I_[model_column][alphabet[A_k[l_count]]].value() + 1;
                }
            }
            std::size_t i = l_count + 1;
            // Inutile de checker la size, la séquence est toujours finie par M.
            // On cherche le prochain état valide (non None)
            while (Pi_k[i] == HMMState::None) {
                i++;
            }
            // Mettre à jour T_u s'il y a un prochain état valable
            if (i < Pi_k.size() && Pi_k[l_count] != HMMState::None) {
                T_[model_column][static_cast<int>(Pi_k[l_count]) * 3 + static_cast<int>(Pi_k[i])] =
                        T_[model_column][static_cast<int>(Pi_k[l_count]) * 3 + static_cast<int>(Pi_k[i])].value() + 1;
            }
            l_count++;
        }
    }
    normalize_matrixes();
    // Ne donne pas de meilleurs résultats, mais pour la forme
    round_matrix(T_);
    round_matrix(e_M_);
    round_matrix(e_I_);
}

// Normalisation pour HMM-build
void HMM::normalize_matrixes() {
    float sum_t;
    float sum_i;
    float sum_m;
    // T : Pour chaque ligne
    for (auto i = 0; i < N_; i++) {
        // Pour chaque état (M, D, I)
        for (auto state = 0; state < 3; state++) {
            sum_t = 0;
            // Calcul de la somme par état (M, D, I)
            for (auto transition = 0; transition < 3; transition++) {
                sum_t += T_[i][3 * state + transition].value();
            }
            // Normalisation par état (M, D, I)
            if (sum_t != 0) {
                for (auto transition = 0; transition < 3; transition++) {
                    T_[i][3 * state + transition].value() =
                            T_[i][3 * state + transition].value() / sum_t;
                }
            }
        }
    }
    // e_M : La somme de chaque ligne vaut 1 (on ignore la première ligne de NAN)
    for (auto line = 1; line < N_; line++) {
        sum_m = 0;
        // Somme
        for (auto column = 0; column < 20; column++) {
            sum_m += e_M_[line][column].value();

        }
        // Normalisation
        for (auto column = 0; column < 20; column++) {
            e_M_[line][column]  = e_M_[line][column].value() / sum_m;

        }
    }
    // e_I : La somme de chaque ligne vaut 1
    for (auto line = 0; line < N_; line++) {
        sum_i = 0;
        // Somme
        for (auto column = 0; column < 20; column++) {
            sum_i += e_I_[line][column].value();

        }
        // Normalisation
        for (auto column = 0; column < 20; column++) {
            e_I_[line][column] = e_I_[line][column].value() / sum_i;

        }
    }
}

// Construction de Pi_k à partir des séquences (HMM-build)
std::vector<HMM::HMMState> HMM::build_Pi_k(const std::vector<char>& sequence) {
    std::vector<HMMState> ret;
    auto column_count = 0;
    for (auto & l : sequence) {
        // Si la colonne est marquée
        if (marked_columns_[column_count]) {
            // Si c'est un '-', alors c'est un état D
            if (l == '-') {
                ret.emplace_back(HMMState::D);
            }
            // Sinon c'est un état M
            else {
                ret.emplace_back(HMMState::M);
            }
        }
        // Si la colonne l n'est pas marquée
        else {
            // Si c'est un '-', alors elle n'a pas d'état
            if (l == '-') {
                ret.emplace_back(HMMState::None);
            }
            // Sinon c'est un état I
            else {
                ret.emplace_back(HMMState::I);
            }
        }
        column_count++;
    }
    // Les états vont finalement vers M, on ajoute donc un état pour modéliser ça
    ret.emplace_back(HMMState::M);
    return ret;
}


// HMM-genseq
void HMM::build_print_genseq() {
    // Les valeurs à écrire en sortie
    std::string sequence;
    std::string states_sequence;
    // On part toujours de l'état M
    HMMState current_state = HMMState::M;
    // Compteurs pour éviter les boucles infinies sur I
    auto k_i = 1;
    bool was_i = false;
    // Initialisation au rang 0 du modèle
    auto chain_index = 0;
    // Pour chaque rang du modèle
    while(chain_index < N_) {
        // Calcul du prochain état : on cherche la probabilité la plus haute de l'ancien état
        current_state = static_cast<HMMState>(index_of_max(T_[chain_index],
                                                           k_i,
                                                           static_cast<int>(current_state) * 3,
                                                           static_cast<int>(current_state) * 3 + 3) % 3);
        // Si l'état est M ou D, incrémenter le rang
        if (current_state == HMMState::M || current_state == HMMState::D) {
            chain_index++;
        }
        // Si nous ne sommes pas au bout du modèle, ajouter l'état actuel à la séquence d'état
        if (chain_index < N_) {
            switch (current_state) {
                case HMMState::M:
                    states_sequence.push_back('M');
                    // Réinitialisation des compteurs de l'état I
                    was_i = false;
                    k_i = 1;
                    break;
                case HMMState::D:
                    states_sequence.push_back('D');
                    // Réinitialisation des compteurs de l'état I
                    was_i = false;
                    k_i = 1;
                    break;
                case HMMState::I:
                    states_sequence.push_back('I');
                    // Si nous étions déjà sur un état I (et donc qu'on boucle sur I), incrémenter le facteur de
                    // puissance à appliquer à la probabilité de transition de I à I
                    if (was_i) {
                        k_i += 1;
                    } else {
                        // Premier état I
                        was_i = true;
                    }
                    break;
                case HMMState::None:
                    // Impossible, pour la forme
                    break;
            }
            // Si c'est un état M, on ajoute à la chaîne le caractère d'émission d'état M le plus probable à ce rang du
            // modèle
            if (current_state == HMMState::M) {
                sequence.push_back(most_probable_char(e_M_[chain_index]));
            } else if (current_state == HMMState::I) { // Similaire pour I
                sequence.push_back(most_probable_char(e_I_[chain_index]));
            } else { // C'est un état D (deletion), on ajoute un '-'
                sequence.push_back('-');
            }
        }
    }
    // Écriture des résultats sur la sortie standard.
    std::cout << sequence << std::endl;
    std::cout << states_sequence << std::endl;
}

// Recherche du caractère le plus probable
char HMM::most_probable_char(std::vector<std::optional<float>> &vector) {
    // Index
    auto max_index = index_of_max(vector);
    // Reverse find
    return find_alphabet_value_of(max_index);
}

// Arrondi
void HMM::round_matrix(std::vector<std::vector<std::optional<float>>> matrix) {
    for (auto & line : matrix) {
        for (auto & element : line) {
            if (element.has_value()) {
                element = round((element.value() * 1000) / 1000);
            }
        }
    }
}

// Setter
void HMM::set_sequences(std::vector<std::vector<char>> sequences) {
    sequences_ = std::move(sequences);
}

// HMM-align
void HMM::viterbi(bool score) {
    // Matrice de score
    std::vector<std::vector<std::optional<float>>> V;
    // Matrice retour
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> B;

    // Modificateurs temporaires sur i et j
    int i_mod;
    int j_mod;
    // Valeurs temporaires pour calculer le maximum et la prochaine valeur de V[i][j], et l'étape retour
    float tmp;
    float v_i_j_value;
    float max_value;
    std::pair<std::size_t, std::size_t> max_coordinates;
    // epsilon
    const float epsilon = 1e-20;

    // Initialisation des matrices sans valeur, sauf la première colonne qui vaut -inf
    for (auto line = 0; line < 3 * N_ + 1; line++) {
        V.emplace_back(sequences_.back().size() + 1, std::optional<float>());
        V.back()[0] = -1 * std::numeric_limits<float>::infinity();
        B.emplace_back(sequences_.back().size() + 1, std::pair<std::size_t, std::size_t>());
    }
    // Première et seconde ligne à -inf
    for (std::size_t column = 0; column < V.back().size(); column ++) {
        V[0][column] = -1 * std::numeric_limits<float>::infinity();
        V[1][column] = -1 * std::numeric_limits<float>::infinity();
    }
    // Sauf V[0][0] = 0
    V[0][0] = 0;
    for(auto i = 2; i < 3 * N_ + 1; i++) {
        for (std::size_t j = 1; j < sequences_.back().size() + 1; j++) {
            // Calcul du prochain état : en fonction de i % 3 : 0 -> M; 1 -> D; 2 -> I (else none pour la forme)
            // Pour chaque état, calcul du modificateur d'index sur i et j ainsi que de la valeur v_i_j propre à chaque
            // état. L'usage des modificateurs permet d'avoir une formule unique dans le calcul du max, car les valeurs
            // à tester sont les mêmes à un facteur constant prêt dans les 3 cas.
            if (i % 3 == 0) /* État M */ {
                i_mod = 1;
                j_mod = 1;
                // Calcul de la "constante" en fonction de e_M
                if (i < 3 * N_) {
                    v_i_j_value = logf(e_M_[i / 3][alphabet[sequences_.back()[j - 1]]].value() + epsilon);
                } else {
                    v_i_j_value = 0;
                    j_mod = 0;
                }

            } else if (i % 3 == 1) /* État D */ {
                i_mod = 2;
                j_mod = 0;
                // Réinitialisation des compteurs
                // La "constante" vaut 0, car c'est un état D
                v_i_j_value = 0;

            } else if (i % 3 == 2) /* État I */ {
                i_mod = 0;
                j_mod = 1;
                // Calcul de la "constante" en fonction de e_I
                v_i_j_value = logf(e_I_[i / 3][alphabet[sequences_.back()[j - 1]]].value() + epsilon);
            }
            // Calculer les 3 valeurs si elles existent (check index)
            // Réinitialisation de la recherche de max
            max_value = -1 * std::numeric_limits<float>::infinity();
            max_coordinates = {};
            // Calcul du maximum des 3 valeurs recherchées. On utilise tmp mod pour itérer sur les 3 valeurs
            for (int tmp_mod = 0; tmp_mod < 3; tmp_mod++) {
                // Calcul de V[][] + log(T[][])
                tmp = V[i - i_mod - tmp_mod][j - j_mod].value() +
                        logf(T_[(i - i_mod - tmp_mod) / 3][(3 * (2 - tmp_mod)) + (i % 3)].value());
                // Comparaison avec le maximum actuel
                if (tmp > max_value) {
                    max_value = tmp;
                    max_coordinates = {i - i_mod - tmp_mod, j - j_mod};
                }
            }
            // Ajouter le maximum au terme d'émission et sauvegarder
            V[i][j] = v_i_j_value + max_value;
            // Étape retour : on retient les coordonnées de la case maximale utilisée pour remplir V[i][j]
            B[i][j] = max_coordinates;
        }
    }
    // Option --score : display le score et quitter
    if (score) {
        std::cout << std::setprecision(3) << std::fixed << V.back().back().value() << std::endl;
        return;
    }
#ifdef DEBUG
    // Debug : print V
    display_matrix(V);
    // Debug : print B
    display_matrix(B);
#endif
    // Étape retour : construction des états et de la séquence alignée : on part de la case en bas à droite de V
    std::pair<std::size_t, std::size_t> current_cell = B.back().back();
    std::string sequence{};
    std::string states_sequence{};
    // Tant que nous n'avons pas atteint un bord (vraisemblablement 0,0)
    while (current_cell.first != 0 && current_cell.second != 0) {
        // Calcul de l'état actuel : la ligne sur laquelle on se trouve % 3 (0=M, 1=D, 2=I)
        // On ajoute l'état actuel au début de la liste d'état (car on reconstruit à l'envers la liste d'états)
        switch (static_cast<HMMState>(current_cell.first % 3)) {
            case HMMState::M:
                states_sequence.insert(0, 1,'M');
                break;
            case HMMState::D:
                states_sequence.insert(0, 1,'D');
                break;
            case HMMState::I:
                states_sequence.insert(0, 1,'I');
                break;
            case HMMState::None:
                break;
        }
        // État suivant
        current_cell = B[current_cell.first][current_cell.second];
    }
    // Reconstruction de la liste d'acides aminés
    auto sequence_index = 0;
    // Pour chaque état
    for (const char & state_count : states_sequence) {
        // Si c'est un M ou I : ajouter la première lettre non déjà écrite
        if (state_count == 'M' || state_count == 'I') {
            sequence.push_back(sequences_.back()[sequence_index++]);
        } else {
            // Sinon, ajouter un '-' pour marquer la deletion
            sequence.push_back('-');
        }
    }
    // output : print les séquences
    std::cout << sequence << std::endl;
    std::cout << states_sequence << std::endl;
}

// Reverse find dans la map (alphabet)
char HMM::find_alphabet_value_of(size_t index) {
    for (auto & mapping : alphabet) {
        if (mapping.second == index) {
            return mapping.first;
        }
    }
    return 0;
}

// Surcharge de display_matrix pour écrire la matrice de paire. Utilisée en débug pour écrire B la matrice retour de
// HMM-align
void HMM::display_matrix(std::vector<std::vector<std::pair<std::size_t, std::size_t>>> matrix) {
    for (auto & line : matrix) {
        for (std::size_t j = 0; j < line.size(); j++) {
                std::cout << "[" << line[j].first << " " << line[j].second << ']';
            if (j != line.size() - 1) {
                std::cout << ',';
            } else {
                std::cout << std::endl;
            }
        }
    }
}