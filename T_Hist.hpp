
class T_Hist
{

public:
    T_Hist();
    T_Hist(std::vector<unsigned char> &);
    SO6 reconstruct();
    T_Hist operator*(T_Hist &);
    bool operator==(T_Hist &);
    bool operator<(T_Hist &);
    bool operator<(const T_Hist &) const;
    friend std::ostream &operator<<(std::ostream &, const T_Hist &);

    std::vector<unsigned char> getHistory() { return hist; };
    void histInsert(unsigned char);

    static void make_T_matrices();
    static SO6 tsv[16];

    struct Node;
    static Node *head;
    static void initHead();
    static void tableInsert(Node *, Node *, unsigned char);
    static void tableDelete(Node *, Node *);
    static SO6 tableLookup(std::vector<unsigned char>);

private:
    std::vector<unsigned char> hist;
};