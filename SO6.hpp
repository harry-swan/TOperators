
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

    static SO6 identity()
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
            Z2 f = first[i];
            Z2 s = second[i];
            if(signA) f.negate();
            if(signB) s.negate();
            if (f != s) return f < s;
        }
        return false;
    }

    static int8_t lexComp(const Z2 first[6], const Z2 second[6])
    {
        for (char i = 0; i < 6; i++)
        {
            if (first[i] < second[i])
                return -1;
            if (second[i] < first[i])
                return 1;
        }
        return 0;
    }

    static int8_t lexComp(Z2 *first, Z2 *second, bool signA, bool signB)
    {
        for (unsigned char i = 0; i < 6; i++)
        {
            Z2 f = first[i];
            Z2 s = second[i];
            if(signA) f.negate();
            if(signB) s.negate();
            if (f < s) return 1;
            if (f > s) return -1;
        }
        return 0;
    }

    static int8_t lexComp(std::vector<Z2> &first, std::vector<Z2> &second, bool signA, bool signB)
    {
        for (unsigned char i = 0; i < 6; i++)
        {
            Z2 f = first[i];
            Z2 s = second[i];
            if(signA) f.negate();
            if(signB) s.negate();
            if (f < s) return 1;
            if (f > s) return -1;
        }
        return 0;
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
            ret[col] = (mat[col][row] < 0);
        }
        return ret;
    }

    static std::vector<int8_t> lexicographic_permutation(SO6 &mat)
    {
        std::vector<bool> signs = column_signs(mat);
        std::vector<int8_t> index(6, 0);
        for (int i = 0; i < 6; i++)
            if(signs[i]) index[i] = -i-1;
            else index[i] = i+1;
        std::sort(index.begin(), index.end(),
                  [&](const int &a, const int &b)
                  {
                      return SO6::lexLess(mat[abs(a)-1], mat[abs(b)-1], signs[abs(a)-1], signs[abs(b)-1]);
                  });
        return index;
    }

    static SO6 const permute_matrix(const std::vector<int> &perms)
    {
        SO6 ret;
        for(int col = 0; col < 6; col++) {
                if (perms[col] > 0) {
                    ret[col][perms[col]-1] = 1;
                }
                else {
                    ret[col][-perms[col]-1] = -1;
                }
        }
        return ret;
    }

    static std::vector<Z2> multiply_only_column(SO6& first, SO6& second, int8_t & col) {
        std::vector<Z2> ret;
        for(int i = 0; i < 6; i++) ret.push_back(Z2(0,0,0));
        Z2 next;
        for(int row = 0; row < 6; row++) {
                for (char k = 0; k < 6; k++) {
                    next = first[k][row] * second[col][k];
                    ret[row] += next; 
                }
        }
        return ret;
    }

    static void print_mat(SO6 & mat) {
        for (int row = 0; row < 6; row++) {
            for (int col = 0; col < 6; col++) {
                std::cout << "(" << static_cast<int>(mat[col][row].val[0]) << static_cast<int>(mat[col][row].val[1]) << static_cast<int>(mat[col][row].val[2]) << ") ";
            }
            std::cout << "\n";
        }
    }

    // static const SO6 permute_matrix(std::vector<int> &perms)
    // {
    //     SO6 ret;
    //     for(int col = 0; col < 6; col++) {
    //             if (perms[col] > 0) {
    //                 ret[col][perms[col]-1] = 1;
    //             }
    //             else {
    //                 ret[col][-perms[col]-1] = -1;
    //             }
    //     }
    //     return ret;
    // }

// private:
    Z2 arr[6][6];
    // std::vector<char> hist;
    // std::string name;
    // Z2 norm;
    // char LDE;
};
