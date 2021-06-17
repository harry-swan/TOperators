class T_Hist 
{

public:
    T_Hist();
    T_Hist(SO6 &);
    T_Hist(std::vector<int8_t> &);
    SO6 reconstruct() ;
    SO6 reconstruct(std::vector<int8_t> &) ;
    SO6 lex_reconstruct() const;
    void set_perm();
    T_Hist operator*(T_Hist &);
    bool operator<(const T_Hist &) const;
    std::vector<int8_t> getHistory() { return hist; };
    static SO6* tsv;

    // static SO6 eval(std::vector<int8_t> &history) {
    //     if(history.size()==0) return tsv[0];
    //     if(history.size()==1) return tsv[history[0]];
    //     std::vector<int8_t> left(history.begin(),history.begin()+history.size()/2);
    //     std::vector<int8_t> right(history.begin()+history.size()/2,history.end());
    //     SO6 ret = eval(left);
    //     SO6 ret2 = eval(right);
    //     ret = ret*ret2;
    //     return ret;
    // }


private:
    std::vector<int8_t> hist;
    std::vector<int> perm;
};