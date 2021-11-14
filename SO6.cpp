#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <bitset>
#include <set>
#include <iterator>
#include <stdint.h>
#include "Z2.hpp"
#include "SO6.hpp"

/**
 * @brief Method to avoid multiple calls to lexLess when we need to lexicographically compare two strings
 * It doesn't seem like this will hit 0, since this only used when sorting, I think...
 * 
 * @param first 
 * @param second 
 * @return char 
 */
// int lexComp (const Z2 first[6], const Z2 second[6]) {
//     for(char i = 0; i < 6 ; i++) {
//         if(first[i] < second[i]) return 1;
//         if(second[i] < first[i]) return -1;
//     }
//     return 0;
// }

// unsigned long long SO6::calls[6] = {0,0,0,0,0,0};

/**
 * Basic constructor. Initializes Zero matrix.
 *
 */
SO6::SO6()
{
    // for(char i=0; i<6; i++){
    //     for(char j=0; j<6; j++)
    //         arr[i][j]=Z2();
    // }
}

/**
 * Constructor that initializes arbitrary matrix with arbitrary name
 * @param a array of Z2 that the SO6 will take as values
 * @param t the object history
 */
SO6::SO6(Z2 a[6][6])
{
    // initializes SO6's entries according to a
    for (unsigned char col = 0; col < 6; col++)
    {
        for (unsigned char row = 0; row < 6; row++)
        {
            arr[col][row] = a[row][col];
        }
    }
}

//Functions to ensure patterns are properly sorted

char SO6::resLexComp(std::vector<char> &first, std::vector<char> &second){
	if(first[0] < second[0]) return 1;
	else if(first[0] > second[0]) return -1;
	else if(first[1] < second[1]) return 1;
	else if(first[1] > second[1]) return -1;
	return 0;
}

bool SO6::resVectLexComp(std::vector<std::vector<char>> &first, std::vector<std::vector<char>> &second){
	char c;
	for(unsigned char i = 0; i<6; i++){
	   c = SO6::resLexComp(first[i], second[i]);
	   if(c==1) return true;
	   else if(c==-1) return false;
    }
    return false;

}

std::vector<char> SO6::resColSort(std::vector<std::vector<std::vector<char>>> &mat){
	std::vector<char> index= {0,1,2,3,4,5};
    std::sort(index.begin(), index.end(),[&](const int &a, const int &b){return SO6::resVectLexComp(mat[a], mat[b]);});
	return index;	
}

std::vector<std::vector<std::vector<char>>> SO6::res_sort(std::vector<std::vector<std::vector<char>>> &mat){
	std::vector<char> permuted = SO6::resColSort(mat);
	std::vector<std::vector<std::vector<char>>> to_return;
	to_return.resize(6);
	std::vector<std::vector<char>> col = mat[0];
	for(char i = 0; i<6; i++){
		to_return[i].resize(6);
		to_return[i] = mat[permuted[i]];
	}
	return to_return;
}



/**
 * Overloads the * operator with matrix multiplication for SO6 objects
 * @param other reference to SO6 to be multiplied with (*this)
 * @return matrix multiplication of (*this) and other
 */
SO6 SO6::operator*(SO6 &other)
{
    //multiplies operators assuming COLUMN,ROW indexing
    SO6 prod;
    Z2 next;

    // Compute product
    for (unsigned char row = 0; row < 6; row++)
    {
        for (unsigned char col = 0; col < 6; col++)
        {
            for (unsigned char k = 0; k < 6; k++)
            {
                // next = arr[row][k]*other[col][k];            // This transpose * other
                next = arr[k][row] * other[col][k]; // This not transpose * other
                // prod(col,row) += next;
                prod[col][row] += next;
            }
        }
    }
    return prod;
}

SO6 SO6::operator*(const SO6 &other) const
{
    //multiplies operators assuming COLUMN,ROW indexing
    SO6 prod;
    Z2 next;

    // Compute product
    for (unsigned char row = 0; row < 6; row++)
    {
        for (unsigned char col = 0; col < 6; col++)
        {
            for (unsigned char k = 0; k < 6; k++)
            {
                // next = arr[row][k]*other[col][k];            // This transpose * other
                next = arr[k][row] * other[col][k]; // This not transpose * other
                // prod(col,row) += next;
                prod[col][row] += next;
            }
        }
    }
    return prod;
}

SO6 SO6::transpose()
{
    return SO6(arr);
}

void SO6::fixSign()
{
    for (unsigned char col = 0; col < 6; col++)
    {
        for (unsigned char row = 0; row < 6; row++)
        {
            if (arr[col][row] < 0)
            {
                while (row < 6)
                {
                    arr[col][row] = -arr[col][row];
                    row++;
                }
            }
            else if (arr[col][row] == 0)
                continue;
            break;
        }
    }
}

// This may be slow
void SO6::lexOrder()
{
    Z2 *myZ2[] = {arr[0], arr[1], arr[2], arr[3], arr[4], arr[5]};
    std::vector<Z2 *> myvector(myZ2, myZ2 + 6);
    std::sort(myvector.begin(), myvector.end(),
              [&](Z2 *a, Z2 *b)
              {
                  return SO6::lexLess(a, b);
              });
    Z2 arr2[6][6];
    for (unsigned char i = 0; i < 6; i++)
    {
        for (unsigned char j = 0; j < 6; j++)
        {
            arr2[i][j] = (myvector.at(i))[j];
        }
    }
    for (unsigned char i = 0; i < 6; i++)
    {
        for (unsigned char j = 0; j < 6; j++)
        {
            arr[i][j] = arr2[i][j];
        }
    }
}

void SO6::reduced_rep()
{
    fixSign();
    lexOrder();
}

bool SO6::operator<(const SO6 &other) const
{
    SO6 first = *this;
    SO6 second = other;

    for (unsigned char col = 0; col < 6; col++)
    {
        switch (lexComp(first[col], second[col]))
        {
        case -1:
            return true;
        case 1:
            return false;
        }
    }
    return false;
}

//std::vector<Z2> getCol(const char &i) {
//    
//}
/** overloads == method to check equality of SO6 matrices
 *  @param other reference to SO6 to be checked against
 *  @return whether or not (*this) and other are equivalent
 */
bool SO6::operator==(SO6 &other)
{
    SO6 first = *this;

    SO6 second = other;
    // SO6 are the same if they have the same triangle
    // TODO: lower right triangle seems super fast, but can try out others
    for (unsigned char col = 5; col != 255; col--)
    {
        for (unsigned char row = 5; row != 255; row--)
        {
            if (first[col][row] != second[col][row])
                return false;
        }
    }
    return true;
}

bool SO6::operator==(const SO6 &other) const
{
    SO6 first = *this;

    SO6 second = other;

    // SO6 are the same if they have the same triangle
    // TODO: lower right triangle seems super fast, but can try out others
    for (unsigned char col = 5; col != 255; col--)
    {
        for (unsigned char row = 5; row != 255; row--)
        {
            if (first[col][row] < second[col][row] || second[col][row] < first[col][row])
                return false;
        }
    }
    return true;
}

/**
 * Overloads << function for SO6.
 * @param os reference to ostream object needed to implement <<
 * @param m reference to SO6 object to be displayed
 * @returns reference ostream with the matrix's display form appended
 */
std::ostream &operator<<(std::ostream &os, const SO6 &m)
{
    std::vector<std::vector<std::vector<char>>> pat = ((SO6)m).pattern();
    os << std::hex << static_cast<int>(((SO6)m).genLDE());
    for(int row = 0; row < 6; row++)
    {
        os << '[';
        for(int col = 0; col < 6; col++)
        {
            os << std::hex << static_cast<int>(pat[col][row][0]) << ' ' << static_cast<int>(pat[col][row][1]) << (col!=5?",":"");
        }
        os << "]";
    }
    return os;
}

char SO6::genLDE()
{
    char LDE = -1;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            if (arr[i][j].getLDE() > LDE)
                LDE = arr[i][j].getLDE();
        }
    }
    return LDE;
}

std::vector<std::vector<std::vector<char>>> SO6::pattern()
{
    char LDE = genLDE();
    std::vector<std::vector<std::vector<char>>> pat;
    pat.resize(6);
    // res.hist = hist;
    for (unsigned char i = 0; i < 6; i++)
    {
        pat[i].resize(6);
        for (unsigned char j = 0; j < 6; j++)
        {
            std::vector<unsigned char> res = arr[i][j].residue(LDE);
            for (char k = 1; k >= 0; k--)
                pat[i][j].emplace_back(res[k]);
        }
    }
    return res_sort(pat);
}

// char lexComp(std::vector<Z2> &first, std::vector<Z2> &second, bool &signA, bool &signB)
//     {
//         for (unsigned char i = 0; i < 6; i++)
//         {
//             Z2 f = first[i];
//             Z2 s = second[i];
//             if(signA) f.negate();
//             if(signB) s.negate();
//             if (f < s) {
//                 SO6::calls.insert(i);
//                 return 1;
//             }
//             if (f > s) {
//                 SO6::calls.insert(i);
//                 return -1;
//             }
//         }
//         SO6::calls.insert(5);
//         return 0;
//     }
