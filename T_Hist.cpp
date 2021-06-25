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

// Nodes of a tree that stores SO6 matrices
// next[i]->so6 = tsv[i+1] * so6
struct T_Hist::Node
{
    SO6 so6;
    Node *next[15];
    Node *prev;
};

SO6 T_Hist::tsv[16] = {SO6::identity(), SO6::tMatrix(0, 1, 0), SO6::tMatrix(0, 2, 1), SO6::tMatrix(0, 3, 2), SO6::tMatrix(0, 4, 3), SO6::tMatrix(0, 5, 4), SO6::tMatrix(1, 2, 5), SO6::tMatrix(1, 3, 6), SO6::tMatrix(1, 4, 7), SO6::tMatrix(1, 5, 8), SO6::tMatrix(2, 3, 9), SO6::tMatrix(2, 4, 10), SO6::tMatrix(2, 5, 11), SO6::tMatrix(3, 4, 12), SO6::tMatrix(3, 5, 13), SO6::tMatrix(4, 5, 14)};
T_Hist::Node *T_Hist::head = NULL;

T_Hist::T_Hist()
{
    hist.clear();
    if(perm.empty()) exit(0);
    SO6 tmp = reconstruct();
    perm = SO6::lexicographic_permutation(tmp);
}

T_Hist::T_Hist(std::vector<unsigned char> &new_hist)
{
    hist = new_hist;
    SO6 tmp = reconstruct();
    perm = SO6::lexicographic_permutation(tmp);
}

void T_Hist::initHead()
{
    T_Hist::head = new Node;
    T_Hist::head->so6 = SO6::identity();
}

// Recursively populates the so6 tree up to depth multiplications
void T_Hist::tableInsert(Node *t, Node *p, unsigned char depth)
{
    if (depth)
    {
        unsigned char i = 0;
        while (i < 15)
        {
            Node *node = new Node;
            t->next[i] = node;
            node->prev = p;
            node->so6 = t->so6 * T_Hist::tsv[++i];
            tableInsert(node, t, depth - 1);
        }
    }
}

// Recursively frees all memory used allocated by the so6 tree
void T_Hist::tableDelete(Node *t, Node *p)
{
    if (t)
    {
        for (unsigned char i = 0; i < 15; i++)
        {
            tableDelete(t->next[i], t);
        }
        delete t;
    }
}

const SO6 T_Hist::tableLookup(std::vector<unsigned char> index)
{
    Node *node = T_Hist::head;
    for (char i : index)
    {
        node = node->next[i - 1];
    }
    return node->so6;
}

/**
 * Reconstructs the corresponding SO6 object by multiplying the left and right SO6 objects
 * found from looking up the left and right history vectors with tableLookup()
 * @param t The SO6 object resulting from multiplication of every tsv element in hist
 */
SO6 T_Hist::reconstruct()
{
    if(hist.size()==0) return SO6::identity();
    std::vector<unsigned char> left(hist.begin(), hist.begin() + hist.size() / 2);
    std::vector<unsigned char> right(hist.begin() + hist.size() / 2, hist.end());
    SO6 ret = tableLookup(left) * tableLookup(right);
    return ret;
}

SO6 T_Hist::reconstruct() const{
    std::vector<unsigned char> left(hist.begin(), hist.begin() + hist.size() / 2);
    std::vector<unsigned char> right(hist.begin() + hist.size() / 2, hist.end());
    SO6 ret = tableLookup(left) * tableLookup(right);
    return ret;
}

std::vector<Z2> T_Hist::reconstruct_col(int & col) const{
    std::vector<unsigned char> left(hist.begin(), hist.begin() + hist.size() / 2);
    std::vector<unsigned char> right(hist.begin() + hist.size() / 2, hist.end());
    SO6 first = tableLookup(left);
    SO6 second = tableLookup(right);
    std::vector<Z2> ret = SO6::multiply_only_column(first,second,col);
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
    const T_Hist first = *this;
    const T_Hist second = other;
    return first < second;
}

// Compares the SO6 objects corresponding to *this and other using the lexicographic ordering
bool T_Hist::operator<(const T_Hist &other) const
{
    // We can maybe speed this up by only reconstructing one column at a time. 
    // This should be possible since we're only over multiplying together two matrices
    // SO6 t = reconstruct();
    // SO6 o = other.reconstruct();
    bool signT,signO;
    int colT, colO;
    for(int col = 0; col<6 ; col++) {
        
        signT = perm[col]<0;
        colT = abs(perm[col])-1;
        std::vector<Z2> tcol = this->reconstruct_col(colT);             // Only need the column we need
        

        signO = other.perm[col]<0;
        colO = abs(other.perm[col])-1;
        std::vector<Z2> ocol = other.reconstruct_col(colO);             // Only need the column we need

        int8_t tmp = SO6::lexComp(tcol,ocol,signT,signO); // For whatever reason using lexLess here gives issues.
        // int8_t tmp = SO6::lexComp(t[colT],o[colO],signT,signO); // For whatever reason using lexLess here gives issues.
        if(tmp != 0) return tmp < 0;
    } 
    return false;
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