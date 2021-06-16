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
 }

T_Hist::T_Hist(SO6 &so6) {
    hist = so6.getHistory();
 }

T_Hist::T_Hist(std::vector<int8_t> & new_hist) {
    hist = new_hist;
}

SO6 T_Hist::reconstruct() {
    SO6 tmp = SO6::identity({});
    std::vector<SO6> to_multiply;
    for(int8_t i : hist) if(!(i==0)) to_multiply.push_back(tsv[i]);                    // Ignore identity         
    std::reverse(to_multiply.begin(),to_multiply.end());                            // Reverse the order to make sure multiplication is handled appropriately
    for(SO6 next : to_multiply) tmp = next*tmp;                                     // Multiply the sequence
    return tmp;
}

T_Hist T_Hist::operator*(T_Hist & other) {
    std::vector<int8_t> history;
    for(int8_t i : hist) if(i!=0) history.push_back(i);
    for(int8_t i : other.hist) if(i!=0) history.push_back(i);
    return T_Hist(history);
}

static SO6 tmp[] = {SO6::identity({}),SO6::tMatrix(0, 1, 1),SO6::tMatrix(0, 2, 2),SO6::tMatrix(0, 3, 3),SO6::tMatrix(0, 4, 4),SO6::tMatrix(0, 5, 5),SO6::tMatrix(1, 2, 6),SO6::tMatrix(1, 3, 7),SO6::tMatrix(1, 4, 8),SO6::tMatrix(1, 5, 9),SO6::tMatrix(2, 3, 10),SO6::tMatrix(2, 4, 11),SO6::tMatrix(2, 5, 12),SO6::tMatrix(3, 4, 13),SO6::tMatrix(3, 5, 14),SO6::tMatrix(4, 5, 15)};
SO6* T_Hist::tsv = tmp;