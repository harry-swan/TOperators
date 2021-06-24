#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <bitset>
#include <iterator>
#include <numeric>
#include "Z2.hpp"
#include "SO6.hpp"
#include "T_Hist.hpp"

SO6 T_Hist::tsv[16] = {SO6::identity(), SO6::tMatrix(0, 1, 0), SO6::tMatrix(0, 2, 1), SO6::tMatrix(0, 3, 2), SO6::tMatrix(0, 4, 3), SO6::tMatrix(0, 5, 4), SO6::tMatrix(1, 2, 5), SO6::tMatrix(1, 3, 6), SO6::tMatrix(1, 4, 7), SO6::tMatrix(1, 5, 8), SO6::tMatrix(2, 3, 9), SO6::tMatrix(2, 4, 10), SO6::tMatrix(2, 5, 11), SO6::tMatrix(3, 4, 12), SO6::tMatrix(3, 5, 13), SO6::tMatrix(4, 5, 14)};
std::vector<SO6> T_Hist::table = {};
std::vector<unsigned long long> T_Hist::offset = {};

T_Hist::T_Hist()
{
    hist.clear();
}

T_Hist::T_Hist(std::vector<unsigned char> &new_hist)
{
    hist = new_hist;
}

void T_Hist::initTable(unsigned char depth)
{
    offset.push_back(1);
    for (unsigned char i = 1; i < depth; i++)
    {
        offset.push_back(offset.back() * 16);
    }
    unsigned long long end = offset.back() * 16;
    for (unsigned long long i = 0; i < end; i++)
    {
        SO6 product = SO6::identity();
        for (unsigned char j = 0; j < depth; j++)
        {
            product = product * tsv[(i / offset[j]) % 16];
        }
        table.push_back(product);
    }
}

SO6 T_Hist::tableLookup(std::vector<unsigned char> index)
{
    unsigned long long product = 0;
    unsigned long long end = index.size();
    for (unsigned long long i = 0; i < end; i++)
    {
        product += index[i] * offset[i];
    }
    return table[product];
}

/**
 * Reconstructs the corresponding SO6 object by multiplying the left and right SO6 objects
 * found from looking up the left and right history vectors with tableLookup()
 * @param t The SO6 object resulting from multiplication of every tsv element in hist
 */
SO6 T_Hist::reconstruct()
{
    std::vector<unsigned char> left(hist.begin(), hist.begin() + hist.size() / 2);
    std::vector<unsigned char> right(hist.begin() + hist.size() / 2, hist.end());
    SO6 ret = tableLookup(left) * tableLookup(right);
    ret.reduced_rep();
    return ret;
}

// Does nothing, will be useful if the logic is needed in multiple functions
// such as * and the constructor
void T_Hist::histInsert(unsigned char h)
{
    // Might want this to be a function
}

// Concatenates the history vectors into one T_Hist object
T_Hist T_Hist::operator*(T_Hist &other)
{
    std::vector<unsigned char> history;
    for (unsigned char i : hist)
        history.push_back(i);
    for (unsigned char i : other.hist)
        history.push_back(i);
    return T_Hist(history);
}

bool T_Hist::operator==(T_Hist &other)
{
    return this->reconstruct() == other.reconstruct();
}

bool T_Hist::operator<(T_Hist &other)
{
    return this->reconstruct() < other.reconstruct();
}

// Compares the SO6 objects corresponding to *this and other using the lexicographic ordering
bool T_Hist::operator<(const T_Hist &other) const
{
    T_Hist t = *this;
    T_Hist o = other;
    return t.reconstruct() < o.reconstruct();
}

// Prints every element of the history vector for T_Hist h
std::ostream &operator<<(std::ostream &os, const T_Hist &h)
{
    for (char i : h.hist)
    {
        os << std::hex << +(i);
    }
    return os;
}