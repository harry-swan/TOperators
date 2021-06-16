class T_Hist 
{

public:
    T_Hist();
    T_Hist(SO6 &);
    T_Hist(std::vector<int8_t> &);
    SO6 reconstruct();
    T_Hist operator*(T_Hist &);

    std::vector<int8_t> getHistory() { return hist; };

    static void make_T_matrices();
    static SO6* tsv;

private:
    std::vector<int8_t> hist;
};