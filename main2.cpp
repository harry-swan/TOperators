/**
 * T Operator Product Generation Main File
 * @file main.cpp
 * @author Swan Klein
 * @author Connor Mooney
 * @author Michael Jarret
 * @author Andrew Glaudell
 * @author Jacob Weston
 * @author Mingzhen Tian
 * @version 6/12/21
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <thread>
#include <future>
#include <fstream>
#include <chrono>
#include <set>
#include <string>
#include <sstream>
#include <functional>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include "Z2.hpp"
#include "SO6.hpp"
#include "T_Hist.hpp"

using namespace std;

const unsigned char numThreads = 50;
unsigned long long operationsPerThread;
unsigned char rem;

const unsigned char tCount = 10;
const Z2 inverse_root2 = Z2::inverse_root2();
const Z2 one = Z2::one();

//Turn this on if you want to read in saved data
const bool tIO = false;
//If tIO true, choose which tCount to begin generating from:
const unsigned char genFrom = tCount;

//Saves every saveInterval iterations
//This also determines parallel block sizes
const unsigned int saveInterval = 50000;

// SO6 identity()
// {
//     SO6 I = SO6::identity();
//     I.lexOrder();
//     return I;
// }

/**
 * Returns the SO6 object corresponding to T[i+1, j+1]
 * @param i the first index to determine the T matrix
 * @param j the second index to determine the T matrix
 * @param matNum the index of the SO6 object in the base vector
 * @return T[i+1,j+1]
 */
const SO6 tMatrix(unsigned char i, unsigned char j, char matNum) { return SO6::tMatrix(i, j, matNum); }

// File reading is no functional at the moment, it would need to be made to work with the T_Hist changes
set<T_Hist> fileRead(unsigned char tc)
{
    ifstream tfile;
    tfile.open(("data/T" + to_string(tc) + ".txt").c_str());
    if (!tfile)
    {
        cout << "File does not exist.\n";
        exit(1);
    }
    set<T_Hist> tset;
    char hist;
    unsigned long long i = 0;
    vector<unsigned char> tmp;
    SO6 m;
    while (tfile.get(hist))
    {
        //Convert hex character to integer
        tmp.push_back((hist >= 'a') ? (hist - 'a' + 10) : (hist - '0'));
        if (++i % tc == 0)
        {
            //Commented out error inducing lines
            //m = tsv.at(tmp.at(tmp.size() - 1));
            for (unsigned char k = tmp.size() - 1; k > -1; k--)
            {
                //m = tbase.at(tmp.at(k)) * m;
            }
            //tset.insert(m);
            tmp.clear();
            tfile.get(hist);
        }
    }
    return tset;
}

// This appends set append to the T file
void writeResults(unsigned char i, unsigned long long save, set<T_Hist> &append)
{
    auto start = chrono::high_resolution_clock::now();
    string fileName = "data/T" + to_string(i + 1) + ".sav";
    fstream write = fstream(fileName, std::ios_base::out);
    write << save;
    write.close();
    fileName = "data/T" + to_string(i + 1) + ".txt";
    write = fstream(fileName, std::ios_base::app);
    for (T_Hist n : append)
        write << n;
    write.close();
    auto end = chrono::high_resolution_clock::now();
    auto ret = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << ">>>Wrote T-Count " << (i + 1) << " to 'data/T" << (i + 1) << ".txt' in " << ret << "ms\n";
}

void threadMult(vector<T_Hist> &threadVector, const unsigned char threadNum, const unsigned long long threadContinue,
                const set<T_Hist> &prior, const set<T_Hist> &current)
{
    // Setting up thread iteration, thread multiplies everything in range (titr, citr) to (tend, cend)
    // by all 15 base T Operators
    unsigned long long start = threadContinue + threadNum * operationsPerThread + min(threadNum, rem);
    if (start >= current.size() * 15)
        return;
    unsigned long long end = min(start + operationsPerThread + (threadNum < rem), (unsigned long long)current.size() * 15);
    set<T_Hist>::iterator citr = current.begin();
    set<T_Hist>::iterator cend = current.begin();
    // Add 1 because we don't care about multiplication by identity tsv[0]
    unsigned char t = 1 + start / current.size();
    advance(citr, start % current.size());
    unsigned char tend = 1 + end / current.size();
    advance(cend, end % current.size());

    T_Hist m, curr;
    while (t < 16)
    {
        while (citr != current.end())
        {
            if (t == tend && citr == cend)
                break;
            curr = *citr;
            vector<unsigned char> vec = {t};
            m = T_Hist(vec) * curr;
            if (prior.find(m) == prior.end())
            {
                threadVector.push_back(m);
            }
            citr++;
        }
        if (t == tend)
            break;
        t++;
        citr = current.begin();
    }
}

int main()
{




    auto tbefore = chrono::high_resolution_clock::now();

    operationsPerThread = saveInterval / numThreads;
    rem = saveInterval % numThreads;
    vector<vector<T_Hist>> threadVectors = {};
    for (unsigned int i = 0; i < numThreads; i++)
    {
        threadVectors.push_back(vector<T_Hist>());
    }

    set<T_Hist> prior({});
    set<T_Hist> current({T_Hist()});
    set<T_Hist> next({});
    set<T_Hist> append({});
    ifstream tfile;
    unsigned char start = 0;

    // fileRead() does not work at present
    /* if (tIO && genFrom > 2)
    {
        prior = fileRead(genFrom - 2, tsv);
        current = fileRead(genFrom - 1, tsv);
        start = genFrom - 1;
    } */

    // Initialize the head of the SO6 tree
    T_Hist::initHead();
    // Generate the lookup table
    std::cout << "\nBeginning Table Generation with Depth " << (+tCount + 1) / 2 << "\n";

    T_Hist::tableInsert(T_Hist::head, NULL, (unsigned char)((+tCount + 1) / 2));

    // /* Here's a thing I'm toying with
    // *  >>>>>>>>>>>>>>>>>>>>>>>>
    // *  >>>>>>>>>>>>>>>>>>>>>>>>
    // *  >>>>>>>>>>>>>>>>>>>>>>>>
    // */ 

    std::cout << "Starting Test 1:\n";
    bool test = true;
    for(int tttt = 0; tttt < 100000; tttt++) {
        std::vector<unsigned char> vec; 
        for(int m = 0; m < 5; m++) {
            unsigned char res = 1 + (rand() %15);
            vec.push_back(res);
        }
        T_Hist tmp = T_Hist(vec);
        SO6 first = tmp.reconstruct();
        SO6 second = first;
        second.reduced_rep();
        std::vector<int8_t> perm = SO6::lexicographic_permutation(first);
        
        for(int row = 0; row<6; row++) {
            for(int col = 0; col <6; col ++) {
                int c = perm[col];
                if(c < 0) {test = -first[-c-1][row]==second[col][row];}
                else {test = first[c-1][row]==second[col][row];}
                if(!test) break;
            }
            if(!test) break;
        }
        if(!test) break;

        for(int row = 0; row<6; row++) {
                for(int col = 0; col <6; col ++) {
                    int c = perm[col];
                    bool sign = c <0;
                    test = !(SO6::lexLess(first[abs(c)-1],second[col],sign,0) || SO6::lexLess(second[col],first[abs(c)-1],0,sign)) ;
                    if(!test) break;
                }
                if(!test) break;
        }    
        if(!test) break;


    }

    
    if(test) std::cout << "Successful.\n";
    else {std::cout << "Failed.\n";}
     
     
    //timing
    // Get every T operator
    for (unsigned char i = start; i < tCount; i++)
    {
        std::cout << "\nBeginning T-Count " << (i + 1) << "\n";

        auto start = chrono::high_resolution_clock::now();

        // Main loop here
        unsigned long long save = 0;
        tfile.open(("data/T" + to_string(i + 1) + "sav").c_str());
        if (!tIO || !tfile)
        {
            remove(("data/T" + to_string(i + 1) + ".txt").c_str());
        }
        else
        {
            next = fileRead(i + 1);
            string str;
            getline(tfile, str);
            stringstream s(str);
            getline(s, str, ' ');
            unsigned char save = stoi(str);
        }
        tfile.close();

        unsigned long long saveEnd = current.size() * 15;
        while (save < saveEnd)
        {
            vector<thread> threads = {};
            for (unsigned char i = 0; i < numThreads; i++)
            {
                threads.emplace_back(thread(threadMult, ref(threadVectors[i]), i, save, ref(prior), ref(current)));
            }
            // Do matrix multiplication in threads
            for (unsigned char i = 0; i < numThreads; i++)
            {
                threads[i].join();
                vector<T_Hist>::iterator itr = threadVectors[i].begin();
                vector<T_Hist>::iterator end = threadVectors[i].end();
                while (itr != end)
                {
                    if (next.insert(*itr).second)
                    {
                        // This is only called if the insert into next succeeded
                        append.insert(*itr);
                    }
                    itr++;
                }
                threadVectors[i].clear();
            }
            save += saveInterval;
            writeResults(i, save, append);
            append.clear();
        }

        // End main loop
        auto end = chrono::high_resolution_clock::now();
        // Begin reporting
        auto ret = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        std::cout << ">>>Found " << next.size() << " new matrices in " << ret << "ms\n";
        prior.clear();
        prior.swap(current); // T++
        current.swap(next);  // T++
    }
    // Free all table memory
    T_Hist::tableDelete(T_Hist::head, NULL);
    chrono::duration<double> timeelapsed = chrono::high_resolution_clock::now() - tbefore;
    std::cout << "\nTotal time elapsed: " << chrono::duration_cast<chrono::milliseconds>(timeelapsed).count() << "ms\n";
    return 0;
}















 












    // std::cout << "Test 3:\n";
    // test = true;
    // for(unsigned char i=1;i<16;i++) {
    //     std::vector<unsigned char> v1 = {i};
    //     T_Hist t1 = T_Hist(v1);
    //     for(unsigned char j=1;j<16;j++) {
    //         std::vector<unsigned char> v2 = {j};
    //         T_Hist t2 = T_Hist(v2);
    //         if(!(t1 < t2 || t2 < t1) && i!=j) {
    //             // This shouldn't ever happen
    //             for(int i : t1.getHistory()) std::cout << i << ",";
    //             for(int j : t2.getHistory()) std::cout << j;
    //             std::cout << ":\n\n";
    //             SO6 s1 = t1.reconstruct();
    //             SO6 s= s1;
    //             for(int row = 0; row < 6; row ++) {
    //                 for(int col = 0; col <6 ; col ++) {
    //                     std::cout << "(" << static_cast<int>(s[col][row].val[0]) << static_cast<int>(s[col][row].val[1]) << ") ";
    //                 }
    //                 std::cout << "\n";
    //             }
    //             std::cout << "\n\n";
    //             SO6 s2 = t2.reconstruct();
    //             s=s2;
    //             for(int row = 0; row < 6; row ++) {
    //                 for(int col = 0; col <6 ; col ++) {
    //                     std::cout << "(" << static_cast<int>(s[col][row].val[0]) << static_cast<int>(s[col][row].val[1]) << ") ";
    //                 }
    //                 std::cout << "\n";
    //             }
                
    //             s=s1;
    //             std::cout << "\n\n";
    //             for(int i : SO6::column_signs(s)) std::cout << i << " ";
    //             std::cout <<"\n";
    //             perm = SO6::lexicographic_permutation(s);
    //             for(int i : perm) std::cout << i << " ";
    //             std::cout << "\n\n";
    //             // s = t2.reconstruct();
    //             for(int row = 0; row < 6; row ++) {
    //                 for(int col = 0; col <6 ; col ++) {
    //                             int signT = perm[col]<0;
    //                             int colT = abs(perm[col])-1;
    //                             if(signT) std::cout << "(" << static_cast<int>(-s[colT][row].val[0]) << static_cast<int>(-s[colT][row].val[1]) << ") ";
    //                             if(!signT) std::cout << "(" << static_cast<int>(s[colT][row].val[0]) << static_cast<int>(s[colT][row].val[1]) << ") ";
    //                 }
    //                 std::cout << "\n";
    //             }                
    //             std::cout << "\n\n";
    //             s.reduced_rep();
    //             for(int row = 0; row < 6; row ++) {
    //                 for(int col = 0; col <6 ; col ++) {
    //                             std::cout << "(" << static_cast<int>(s[col][row].val[0]) << static_cast<int>(s[col][row].val[1]) << ") ";                    
    //                 }
    //                 std::cout << "\n";
    //             }                                
                
    //             std::cout << "\n\n\n";

    //             s=s2;
    //             std::cout << "\n\n";
    //             for(int i : SO6::column_signs(s)) std::cout << i << " ";
    //             std::cout <<"\n";
    //             perm = SO6::lexicographic_permutation(s);
    //             for(int i : perm) std::cout << i << " ";
    //             std::cout << "\n\n";
    //             // s = t2.reconstruct();
    //             for(int row = 0; row < 6; row ++) {
    //                 for(int col = 0; col <6 ; col ++) {
    //                             int signT = perm[col]<0;
    //                             int colT = abs(perm[col])-1;
    //                             if(signT) std::cout << "(" << static_cast<int>(-s[colT][row].val[0]) << static_cast<int>(-s[colT][row].val[1]) << ") ";
    //                             if(!signT) std::cout << "(" << static_cast<int>(s[colT][row].val[0]) << static_cast<int>(s[colT][row].val[1]) << ") ";
    //                 }
    //                 std::cout << "\n";
    //             }                
    //             std::cout << "\n\n";
    //             s.reduced_rep();
    //             for(int row = 0; row < 6; row ++) {
    //                 for(int col = 0; col <6 ; col ++) {
    //                             std::cout << "(" << static_cast<int>(s[col][row].val[0]) << static_cast<int>(s[col][row].val[1]) << ") ";                    
    //                 }
    //                 std::cout << "\n";
    //             }                                
                
    //             std::cout << "\n";
    //             test = false;
    //             break;

    //         }

    //         if(!(t1 < t2 || t2 < t1) && i==j) {
    //             for(int i : t1.getHistory()) std::cout << i << " ";
    //             for(int j : t2.getHistory()) std::cout << j << " ";
    //             std::cout << ":\n\n";
    //             SO6 s1 = t1.reconstruct();
    //             SO6 original = s1;
    //             SO6 s2;     // Lexsorted s1 by perm
    //             for(int row = 0; row < 6; row ++) {
    //                 for(int col = 0; col <6 ; col ++) {
    //                     int colT = abs(perm[col])-1;
    //                     s2[colT][row] = s1[col][row];
    //                     if(perm[col] < 0) {
    //                         std::cout << "colT is" << colT << perm[col] << "\n";
    //                         s2[col][row].negate();
    //                     }
    //                 }
    //                 std::cout << "\n";
    //             }                
    //             s1.reduced_rep(); // Lexsort s1
    //             if(!(s1 == s2)) {
    //                 std::cout << "problem!\n";
    //                 for(int i : perm) std::cout << i << " ";
    //                 std::cout << "\n\n";
    //                 for(int row = 0; row < 6; row ++) {
    //                     for(int col = 0; col <6 ; col ++) {
    //                         std::cout << "(" << static_cast<int>(original[col][row].val[0]) << static_cast<int>(original[col][row].val[1]) << ") ";
    //                     }
    //                     std::cout << "\n";
    //                 }
    //                 std::cout << "\n\n";
    //                 for(int row = 0; row < 6; row ++) {
    //                     for(int col = 0; col <6 ; col ++) {
    //                         std::cout << "(" << static_cast<int>(s1[col][row].val[0]) << static_cast<int>(s1[col][row].val[1]) << ") ";
    //                     }
    //                     std::cout << "\n";
    //                 }
    //                 std::cout << "\n\n";
    //                 for(int row = 0; row < 6; row ++) {
    //                     for(int col = 0; col <6 ; col ++) {
    //                         std::cout << "(" << static_cast<int>(s2[col][row].val[0]) << static_cast<int>(s2[col][row].val[1]) << ") ";
    //                     }
    //                     std::cout << "\n";
    //                 }
    //                 test = false;
    //                 break;
    //             }
    //         }
    //         if(!test) break;
    //     }
    //     if(!test) break;
    //     // std::cout << "\n";
    // }
    // if(test) std::cout << "Successful.\n";
    // else std::cout << "Failed.\n";
    // exit(0);