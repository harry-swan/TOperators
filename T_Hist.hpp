
class T_Hist
{

public:
    T_Hist();
    T_Hist(std::vector<int8_t> &);
    SO6 reconstruct();
    T_Hist operator*(T_Hist &);
    bool operator==(T_Hist &);
    bool operator<(T_Hist &);
    bool operator<(const T_Hist &) const;
    friend std::ostream &operator<<(std::ostream &, const T_Hist &);

    std::vector<int8_t> getHistory() { return hist; };

    static void make_T_matrices();
    static SO6 tsv[15];

    struct Node;
    static Node* head;
    static void initHead();
    static void tableInsert(Node *, Node *, int8_t);
    static void tableDelete(Node *, Node *);
    static SO6 tableLookup(std::vector<int8_t>);

private:
    std::vector<int8_t> hist;
};