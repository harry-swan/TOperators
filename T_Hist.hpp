
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

    std::vector<unsigned char> getHistory() { std::vector<unsigned char> ret = left_hist;
                    ret.insert(ret.end(), right_hist.begin(), right_hist.end()); return ret; };
    void histInsert(unsigned char); // Currently unused

    static void make_T_matrices();
    static SO6 tsv[16];

    struct Node;
    static Node *head; // Head of the tree for looking up SO6 objects for reconstruct()
    static void initHead(); // Initalizes the head, called at the start of main()
    static void tableInsert(Node *, Node *, unsigned char); // Populates the so6 tree
    static void tableDelete(Node *, Node *); // Frees all memory allocated in the so6 tree
    static SO6* tableLookup(std::vector<unsigned char> &); // Get the matrix corresponding to a history vector

private:
    std::vector<unsigned char> left_hist;
    std::vector<unsigned char> right_hist;
};