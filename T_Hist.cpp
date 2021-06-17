#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <bitset>
#include <iterator>
#include <stdint.h>
#include "Z2.hpp"
#include "SO6.hpp"
#include "T_Hist.hpp"

T_Hist::T_Hist() {
    hist.clear();
    perm[0] = 7;
 }

T_Hist::T_Hist(SO6 &so6) {
    hist = so6.getHistory();
    perm = SO6::lexicographic_permutation(so6);
 }

T_Hist::T_Hist(std::vector<int8_t> & new_hist) {
    hist = new_hist;
}

SO6 T_Hist::reconstruct() {
    // Recursively multiply, make me better by not declaring new vectors
    std::vector<int8_t> left(hist.begin(),hist.begin()+hist.size()/2);
    std::vector<int8_t> right(hist.begin()+hist.size()/2,hist.end());
    SO6 ret1 = reconstruct(left);
    SO6 ret2 = reconstruct(right);
    return ret1*ret2;
}

// Reconstruct from a given history state. This can probably be made static and should also return a known SO6 if history corresponds to previously evaluated SO6
SO6 T_Hist::reconstruct(std::vector<int8_t> &history) {
    if(history.size()==0) return tsv[0];
    if(history.size()==1) return tsv[history[0]];
    std::vector<int8_t> left(history.begin(),history.begin()+history.size()/2);
    std::vector<int8_t> right(history.begin()+history.size()/2,history.end());
    SO6 ret = reconstruct(left);
    SO6 ret2 = reconstruct(right);
    ret = ret*ret2;
    return ret;
}

void T_Hist::set_perm() {
    SO6 mat = reconstruct();
    perm = SO6::lexicographic_permutation(mat);
 }

SO6 T_Hist::lex_reconstruct() const{
    try {
        if(perm[0]==7) throw (0);
        SO6 mat = reconstruct();
        SO6 ret;
        for(int i = 0; i<6; i++) {
            for(int j = 0; j<6; j++) {
                if(perm[i] > 0) ret[i][j] = mat[perm[i]-1][j];                          // perm is a permutation matrix with entries +/- (1,2,3,4,5,6), need to downshift. maybe a simpler method exists
                else ret[i][j] = -mat[-perm[i]-1][j];                                   // perm is a permutation matrix with entries +/- (1,2,3,4,5,6), need to downshift
            }
        }
        return ret;
    } catch (int n) {
        std::cout << "Warning: comparison between lex-unsorted SO6 objects, sorting, but will be faster if corrected.\n";
        T_Hist th = *this;
        th.set_perm();
        return(th.lex_reconstruct());
    }
}

T_Hist T_Hist::operator*(T_Hist & other) {
    std::vector<int8_t> history;
    for(int8_t i : hist) if(i!=0) history.push_back(i);
    for(int8_t i : other.hist) if(i!=0) history.push_back(i);
    return T_Hist(history);
}

bool T_Hist::operator<(const T_Hist & other) const{
    return (*this).lex_reconstruct() < other.lex_reconstruct(); 
}

static SO6 tmp[] = {SO6::identity({}),SO6::tMatrix(0, 1, 1),SO6::tMatrix(0, 2, 2),SO6::tMatrix(0, 3, 3),SO6::tMatrix(0, 4, 4),SO6::tMatrix(0, 5, 5),SO6::tMatrix(1, 2, 6),SO6::tMatrix(1, 3, 7),SO6::tMatrix(1, 4, 8),SO6::tMatrix(1, 5, 9),SO6::tMatrix(2, 3, 10),SO6::tMatrix(2, 4, 11),SO6::tMatrix(2, 5, 12),SO6::tMatrix(3, 4, 13),SO6::tMatrix(3, 5, 14),SO6::tMatrix(4, 5, 15)};
SO6* T_Hist::tsv = tmp;