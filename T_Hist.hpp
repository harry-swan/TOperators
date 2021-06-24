
class T_Hist
{

public:
    T_Hist();
    T_Hist(std::vector<unsigned char> &);
    SO6 reconstruct(); // Gets an SO6 object by multiplying together the right and left hist
    T_Hist operator*(T_Hist &); // Multiplication is history concatenation
    bool operator==(T_Hist &);
    bool operator<(T_Hist &);
    bool operator<(const T_Hist &) const;
    friend std::ostream &operator<<(std::ostream &, const T_Hist &);

    std::vector<unsigned char> getHistory() { return hist; };
    void histInsert(unsigned char); // Currently unused

    static void make_T_matrices();
    static SO6 tsv[16];

    static std::vector<SO6> table;
    static std::vector<unsigned long long> offset;
    static void initTable(unsigned char);
    static SO6 tableLookup(std::vector<unsigned char>); // Get the matrix corresponding to a history vector

private:
    std::vector<unsigned char> hist;
};