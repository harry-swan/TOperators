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

const unsigned char numThreads = 4;
unsigned long long operationsPerThread;
unsigned char rem;

const unsigned char tCount = 5;
const Z2 inverse_root2 = Z2::inverse_root2();
const Z2 one = Z2::one();

//Turn this on if you want to read in saved data
const bool tIO = false;
//If tIO true, choose which tCount to begin generating from:
const unsigned char genFrom = tCount;

//Saves every saveInterval iterations
//This also determines parallel block sizes
const unsigned int saveInterval = 50000;

SO6 identity()
{
    SO6 I = SO6::identity();
    I.lexOrder();
    return I;
}

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

// This appends next to the T file
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
    // Setting up thread iteration, thread goes from (titr, citr) to (tend, cend)
    unsigned long long start = threadContinue + threadNum * operationsPerThread + min(threadNum, rem);
    if (start >= current.size() * 15)
        return;
    unsigned long long end = min(start + operationsPerThread + (threadNum < rem), (unsigned long long)current.size() * 15);
    set<T_Hist>::iterator citr = current.begin();
    set<T_Hist>::iterator cend = current.begin();
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
    //timing
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

    /* if (tIO && genFrom > 2)
    {
        prior = fileRead(genFrom - 2, tsv);
        current = fileRead(genFrom - 1, tsv);
        start = genFrom - 1;
    } */

    T_Hist::initHead();
    // Generate the lookup table
    std::cout << "\nBeginning Table Generation with Depth " << (+tCount + 1) / 2 << "\n";

    T_Hist::tableInsert(T_Hist::head, NULL, (unsigned char)((+tCount + 1) / 2));

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
    T_Hist::tableDelete(T_Hist::head, NULL);
    chrono::duration<double> timeelapsed = chrono::high_resolution_clock::now() - tbefore;
    std::cout << "\nTotal time elapsed: " << chrono::duration_cast<chrono::milliseconds>(timeelapsed).count() << "ms\n";
    return 0;
}
