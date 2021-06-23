
class SO6
{
public:
    SO6();
    // SO6(std::string); //initializes zero matrix
    // SO6(Z2[6][6], std::string); //initializes matrix according to a 6x6 array of Z2
    SO6(Z2[6][6]); //initializes matrix according to a 6x6 array of Z2
    SO6 operator*(SO6 &);                 //mutliplication
    SO6 operator*(const SO6 &) const;
    void fixSign();
    void lexOrder();
    void reduced_rep();
    inline Z2 &operator()(unsigned char col, unsigned char row) { return arr[col][row]; } //returns the (i,j)th entry
    // bool operator%(SO6&);
    bool operator<(const SO6 &) const;
    const Z2 &operator()(unsigned char i, unsigned char j) const { return arr[i][j]; } //returns the (i,j)th entry but for const
    bool operator==(SO6 &);                                              //checking equality up to signed permutation
    Z2 *operator[](const unsigned char i) { return arr[i]; }                    // Return the array element needed.
    const Z2 *operator[](const unsigned char i) const { return arr[i]; }        // Return the array element needed.
    // inline std::string getName(){return(name);} //getter for Name
    // inline void setName(std::string newName){name = newName;}
    char getLDE() { return (0); } //getter for LDE
    SO6 transpose();
    // void genOneNorm(); //Returns 1-norm of the first row
    // bool isPerm(SO6 s); //Returns true if and only if s is similar to this object
    // inline float normFloat() { return norm.toFloat(); }
    char genLDE();                                              // generates LDE, called after multiplication and constructor
    friend std::ostream &operator<<(std::ostream &, const SO6 &); //display
    bool operator==(const SO6 &other) const;
    SO6 residue(); // returns the residue matrix using Z2::pattern()
    /* inline char getLast()
    {
        if (hist.size() != 0)
            return hist[0];
        return -1;
    }; */
    //std::vector<char> getHistory() { return hist; };

    static const SO6 tMatrix(unsigned char i, unsigned char j, unsigned char matNum)
    {
        SO6 t = SO6::identity(); // Initialize to the identity matrix
        t[i][i] = Z2::inverse_root2();   // Change the i,j cycle to appropriate 1/sqrt(2)
        t[j][j] = Z2::inverse_root2();
        t[i][j] = Z2::inverse_root2();
        if (abs(i - j) != 1)
            t[i][j].negate();
        t[j][i] = -t[i][j];
        if (t == SO6::identity())
        {
            std::cout << static_cast<int16_t>(i) << " " << static_cast<int16_t>(j) << " " << static_cast<int16_t>(matNum) << "\n";
        }
        return t;
    };

    static const SO6 identity()
    {
        SO6 t;
        for (unsigned char k = 0; k < 6; k++)
            t[k][k] = 1;
        return t;
    }

    static bool *signs(SO6 &mat, bool *ret)
    {
        for (unsigned char col = 0; col < 6; col++)
        {
            unsigned char row = 0;
            while (mat[col][row] == 0)
                row++;
            ret[col] = (mat[col][row] < 0);
        }
        return ret;
    }

    // static char* SO6::lexOrder(SO6 &mat)
    // {
    //     char ret[6];
    //     Z2 *myZ2[] = {mat[0], mat[1], mat[2], mat[3], mat[4], mat[5]};
    //     std::vector<Z2 *> sorted(myZ2, myZ2 + 6);
    //     std::sort(sorted.begin(), sorted.end(), SO6::lexLess);
    //     for(char i = 0; i<6; i++) {
    //         if(lexLess(mat[i],sorted.at(2))) {
    //             if(lexLess(mat[i],sorted.at(1))) {
    //                 ret[i] = 0;
    //                 break;
    //             }
    //             ret[i] = 1;
    //             break;
    //         }
    //         else if (lexLess(sorted.at(3),mat[i])) {
    //             if(lexLess(sorted.at(4),mat[i])) {
    //                 ret[i] = 5;
    //                 break;
    //             ret[i] = 1;
    //             break;
    //         }
    //     }
    //     return ret;
    // }

    /**
     * Method to compare two Z2 arrays of length 6 lexicographically
     * @param first array of Z2 of length 6
     * @param second array of Z2 of length 6
     * @return -1 if first < second, 0 if equal, 1 if first > second
     */
    static const bool lexLess(Z2 *first, Z2 *second)
    {
        for (unsigned char i = 0; i < 6; i++)
        {
            if (first[i] != second[i])
                return first[i] < second[i];
        }
        return false;
    }

    static bool lexLess(Z2 *first, Z2 *second, bool signA, bool signB)
    {
        for (unsigned char i = 0; i < 6; i++)
        {
            Z2 f = Z2((1 - 2 * signA), 0, 0);
            Z2 s = Z2((1 - 2 * signB), 0, 0);
            f = f * first[i];
            s = s * second[i];
            if (f != s)
                return f < s;
        }
        return false;
    }

    static std::vector<bool> column_signs(SO6 &mat)
    {
        std::vector<bool> ret(6, 0);
        for (unsigned char col = 0; col < 6; col++)
        {
            unsigned char row = 0;
            while (mat[col][row] == 0)
            {
                row++;
            }
            bool tmp = (mat[col][row] < 0);
            ret[col] = (mat[col][row] < 0);
        }
        return ret;
    }

    static std::vector<int> lexicographic_permutation(SO6 &mat)
    {
        std::vector<bool> signs = column_signs(mat);
        std::vector<int> index(6, 0);
        for (int i = 0; i < 6; i++)
            index[i] = i;
        std::sort(index.begin(), index.end(),
                  [&](const int &a, const int &b)
                  {
                      return SO6::lexLess(mat[a], mat[b], signs[a], signs[b]);
                  });
        for (int i = 0; i < 6; i++)
        {
            if (signs[index[i]])
            {
                index[i] = -(index[i] + 1); // Need negation and cannot negate 0. Can probably be handled without 1 indexing
            }
        }
        return index;
    }

private:
    Z2 arr[6][6];
    // std::vector<char> hist;
    // std::string name;
    // Z2 norm;
    // char LDE;
};
