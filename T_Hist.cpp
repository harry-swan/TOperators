#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <bitset>
#include <iterator>
#include <utility>
#include <stdint.h>
#include "Z2.hpp"
#include "SO6.hpp"
#include "T_Hist.hpp"

// Nodes of a tree that stores SO6 matrices
// next[i]->so6 = tsv[i+1] * so6
struct T_Hist::Node
{
    SO6 *so6;
    Node *next[15];
    Node *prev;
};

SO6 T_Hist::tsv[16] = {SO6::identity(), SO6::tMatrix(0, 1, 0), SO6::tMatrix(0, 2, 1), SO6::tMatrix(0, 3, 2), SO6::tMatrix(0, 4, 3), SO6::tMatrix(0, 5, 4), SO6::tMatrix(1, 2, 5), SO6::tMatrix(1, 3, 6), SO6::tMatrix(1, 4, 7), SO6::tMatrix(1, 5, 8), SO6::tMatrix(2, 3, 9), SO6::tMatrix(2, 4, 10), SO6::tMatrix(2, 5, 11), SO6::tMatrix(3, 4, 12), SO6::tMatrix(3, 5, 13), SO6::tMatrix(4, 5, 14)};
T_Hist::Node *T_Hist::head = NULL;
SO6 *T_Hist::curr = NULL;
T_Hist *T_Hist::curr_history = NULL;

T_Hist::T_Hist()
{
    hist = {0};                         // Initialize to identity
    perm = {6,5,4,3,2,1};               // Identity would get reversed in order
}

T_Hist::T_Hist(unsigned char s)
{
    hist.resize(s, 0);
    perm = {6,5,4,3,2,1};               // I assume here we're making a history of size s all initialized to the identity. 
}

T_Hist::T_Hist(std::vector<unsigned char> &new_hist)
{
    unsigned char s = new_hist.size();
    hist.resize(s%2 + s/2, 0);
    unsigned char i = 0;
    for (unsigned char h : new_hist)
        histInsert(h, i++);

    // I don't know if we do this any longer or, if we do, if it even matters. -Michael
    // We do need this, but only for fileRead(). -Swan
    SO6 tmp = reconstruct();
    perm = SO6::lexicographic_permutation(tmp);
}

void T_Hist::initHead()
{
    T_Hist::head = new Node;
    T_Hist::head->so6 = new SO6(SO6::identity());
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
            node->so6 = new SO6(*(t->so6) * T_Hist::tsv[++i]);
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
        delete t->so6;
        delete t;
    }
}

SO6* T_Hist::tableLookup(std::vector<unsigned char> &index)
{
    Node *node = T_Hist::head;
    if (index.size() == 0)
        return node->so6; // Returns a pointer to the identity
    unsigned char end = 2*index.size() - (index.back() >> 4 == 0);
    // Skip the first 4 bits if they are 0
    unsigned char i = ((index[0] & 15) == 0);
    for (i; i < end; i++)
    {
        //std::cout << +((index[i/2] >> (4*(i%2))) & 15) << "\n";
        node = node->next[((index[i/2] >> (4*(i%2))) & 15) - 1];
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
    unsigned char splitMiddleBit = (hist.size() % 2) && ((hist.back() & 240) != 0);
    std::vector<unsigned char> left(hist.begin(), hist.begin() + hist.size() / 2 + splitMiddleBit);
    std::vector<unsigned char> right(hist.begin() + hist.size() / 2, hist.end());
    if(splitMiddleBit)
    {
        left.back() = left.back() & 15;
        right[0] = right[0] & 240;
    }
    return *tableLookup(left) * *tableLookup(right);
}

std::vector<Z2> T_Hist::reconstruct_col(char & col) const{
    bool splitMiddleBit = (hist.size() % 2) && ((hist.back() & 240) != 0);
    std::vector<unsigned char> left(hist.begin(), hist.begin() + hist.size() / 2 + splitMiddleBit);
    std::vector<unsigned char> right(hist.begin() + hist.size() / 2, hist.end());
    if(splitMiddleBit)
    {
        left.back() = left.back() & 15;
        right[0] = right[0] & 240;
    }
    SO6* first = tableLookup(left);
    SO6* second = tableLookup(right);
    std::vector<Z2> ret = SO6::multiply_only_column(*first,*second,col);
    return ret;
}

void T_Hist::histInsert(unsigned char h, unsigned char i)
{
    hist[i/2] |= h << 4*(i%2);
}

// Concatenates the history vectors into one T_Hist object

T_Hist T_Hist::operator*(T_Hist &other)
{
    // The following lines are commented out since, at this point, no hist.size() should ever be 0. This may simplify some of the logic below and reduce operations elsewhere

    // unsigned char extra = 0;
    // if (hist.size() && other.hist.size())
    //     extra = hist.back() < 16 && other.hist.back() < 16;
    
    unsigned char extra = hist.back() < 16 && other.hist.back() < 16;           

    unsigned char idx = 0;
    T_Hist history(hist.size() + other.hist.size() - extra);

    // if (hist.size())
    // {
    unsigned char end = 2*hist.size() - (hist.back() < 16);
    for (unsigned char i = (hist.back() & 15 == 0); i < end; i++)
        history.histInsert(((hist[i/2] >> (4*(i%2))) & 15), idx++);
    // }
    // if (other.hist.size())
    // {
    end = 2*other.hist.size() - (other.hist.back() < 16);
    for (unsigned char i = (other.hist.back() & 15 == 0); i < end; i++)
        history.histInsert(((other.hist[i/2] >> (4*(i%2))) & 15), idx++);
    // }
    history.set_lex_perm();
    return history;
}

void T_Hist::set_lex_perm() {
    SO6 tmp = reconstruct();
    perm = SO6::lexicographic_permutation(tmp);
}

// Compares the SO6 objects corresponding to *this and other using the lexicographic ordering
bool T_Hist::operator<(const T_Hist &other) const
{
    bool s0,s1;                                // These booleans flag the signs from the permutation list
    char i0,i1;                                // Some integers the columns of interest as flagged by the permutation list
    std::vector<Z2> col0(6), col1(6);           

    if(*this == *T_Hist::curr_history || &other == T_Hist::curr_history) {
        s0 = false;
        SO6 &first = *T_Hist::curr;
        for(int lex_index = 0; lex_index<6; lex_index++) {
            
            for(int row = 0; row < 6; row++) col0[row]=first[lex_index][row];           // May be possible to make faster by doing away with this constructor and passing the array to lexComp below somehow

            if(*this == *T_Hist::curr_history) {
                s1 = other.perm[lex_index]<0;
                i1 = abs(other.perm[lex_index])-1;
                col1 = other.reconstruct_col(i1);
                int8_t ret = SO6::lexComp(col0,col1,s0,s1);                             // For whatever reason using lexLess here gives issues. 
                if(ret != 0) return ret < 0;
            } else {
                s1 = perm[lex_index]<0;
                i1 = abs(perm[lex_index])-1;
                col1 = reconstruct_col(i1);
                int8_t ret = SO6::lexComp(col0,col1,s0,s1);                             // For whatever reason using lexLess here gives issues. 
                if(ret != 0) return ret > 0;
            }
        }
        return false;
    } 

    // If we make it this far, we have to rebuild both matrices. Not sure when/why this happens. Probably called somewhere in the set insert.
    for(int lex_index = 0; lex_index<6 ; lex_index++) {
        s0 = perm[lex_index]<0;
        i0 = abs(perm[lex_index])-1;
        col0 = reconstruct_col(i0);    

        s1 = other.perm[lex_index]<0;
        i1 = abs(other.perm[lex_index])-1;
        col1 = other.reconstruct_col(i1);             

        int8_t tmp = SO6::lexComp(col0,col1,s0,s1); // For whatever reason using lexLess here gives issues. 
        if(tmp != 0) return tmp < 0;
    } 
    return false;
}

// Prints every element of the history vector for T_Hist h
std::ostream &operator<<(std::ostream &os, const T_Hist &h)
{
    if(h.hist.size()==0) return os;                  // This is needed, but I don't know why.
    unsigned char end = 2*h.hist.size() - ((h.hist.back() & 240) == 0);
    for (unsigned char i = 0; i < end; i++)
    {
        os << std::hex << ((h.hist[i/2] >> (4*(i%2))) & 15);
    }
    return os;
}

bool T_Hist::operator==(T_Hist &other)
{
    return &other == this;
}

bool T_Hist::operator==(T_Hist &other) const
{
    return &other == this;
}

// bool T_Hist::operator<(T_Hist &other)
// {
//     const T_Hist first = *this;
//     const T_Hist second = other;
//     return first < second;
// }