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

const int8_t numThreads = 4;
long operationsPerThread;
int8_t rem;

const int8_t tCount = 5;
const Z2 inverse_root2 = Z2::inverse_root2();
const Z2 one = Z2::one();

//Turn this on if you want to read in saved data
const bool tIO = false;
//If tIO true, choose which tCount to begin generating from:
const int8_t genFrom = tCount;

//Saves every saveInterval iterations
//This also determines parallel block sizes
const int saveInterval = 50000;

SO6 identity() 
{
    SO6 I = SO6::identity({0});
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
const SO6 tMatrix(int8_t i, int8_t j, int8_t matNum) {return SO6::tMatrix(i,j,matNum);}

set<SO6> fileRead(int8_t tc, vector<SO6> tbase)
{
    ifstream tfile;
    tfile.open(("data/T" + to_string(tc) + ".txt").c_str());
    if (!tfile)
    {
        cout << "File does not exist.\n";
        exit(1);
    }
    set<SO6> tset;
    char hist;
    long i = 0;
    vector<int8_t> tmp;
    SO6 m;
    while (tfile.get(hist))
    {
        //Convert hex character to integer
        tmp.push_back((hist >= 'a') ? (hist - 'a' + 10) : (hist - '0'));
        if (++i % tc == 0)
        {
            m = tbase.at(tmp.at(tmp.size() - 1));
            for (int8_t k = tmp.size() - 2; k > -1; k--)
            {
                m = tbase.at(tmp.at(k)) * m;
            }
            tset.insert(m);
            tmp.clear();
            tfile.get(hist);
        }
    }
    return tset;
}

// This appends next to the T file
void writeResults(int8_t i, long save, set<SO6> &append)
{
    auto start = chrono::high_resolution_clock::now();
    string fileName = "data/T" + to_string(i + 1) + ".sav";
    fstream write = fstream(fileName, std::ios_base::out);
    write << save;
    write.close();
    fileName = "data/T" + to_string(i + 1) + ".txt";
    write = fstream(fileName, std::ios_base::app);
    for (SO6 n : append)
        write << n;
    write.close();
    auto end = chrono::high_resolution_clock::now();
    auto ret = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << ">>>Wrote T-Count " << (i + 1) << " to 'data/T" << (i + 1) << ".txt' in " << ret << "ms\n";
}

void threadMult(vector<SO6> &threadVector, const int8_t threadNum, const long threadContinue,
                vector<SO6> &tsv, const set<SO6> &prior, const set<SO6> &current)
{
    // Setting up thread iteration, thread goes from (titr, citr) to (tend, cend)
    long start = threadContinue + threadNum * operationsPerThread + min(threadNum, rem);
    if (start >= current.size() * 15)
        return;
    long end = min(start + operationsPerThread + (threadNum < rem), (long)current.size() * 15);
    vector<SO6>::iterator titr = tsv.begin();
    set<SO6>::iterator citr = current.begin();
    vector<SO6>::iterator tend = tsv.begin();
    set<SO6>::iterator cend = current.begin();
    advance(titr, 1 + start / current.size());
    advance(citr, start % current.size());
    advance(tend, 1 + end / current.size());
    advance(cend, end % current.size());

    SO6 t, m, curr;
    while (titr != tsv.end())
    {
        t = *titr;
        while (citr != current.end())
        {
            if (titr == tend && citr == cend)
                break;
            curr = *citr;
            m = t * curr;
            if (prior.find(m) == prior.end())
            {
                threadVector.push_back(m);
            }
            citr++;
        }
        if (titr == tend)
            break;
        titr++;
        citr = current.begin();
    }
}

int main()
{
    operationsPerThread = saveInterval / numThreads;
    rem = saveInterval % numThreads;
    vector<vector<SO6>> threadVectors = {};
    for (int i = 0; i < 15; i++)
    {
        threadVectors.push_back(vector<SO6>());
    }

    //timing
    auto tbefore = chrono::high_resolution_clock::now();

    set<SO6> prior;
    set<SO6> current({identity()});
    set<SO6> next;
    set<SO6> append;
    ifstream tfile;
    int8_t start = 0;

    vector<SO6> tsv(T_Hist::tsv,T_Hist::tsv+16); //t count 1 matrices
    for(int i = 0; i< 16; i++) if(!(tsv[i] == T_Hist::tsv[i])) exit(0); //Failsafe

    if (tIO && genFrom > 2)
    {
        prior = fileRead(genFrom - 2, tsv);
        current = fileRead(genFrom - 1, tsv);
        start = genFrom - 1;
    }

    for (int8_t i = start; i < tCount; i++)
    {
        std::cout << "\nBeginning T-Count " << (i + 1) << "\n";

        auto start = chrono::high_resolution_clock::now();

        // Main loop here
        long save = 0;
        tfile.open(("data/T" + to_string(i + 1) + "sav").c_str());
        if (!tIO || !tfile)
        {
            remove(("data/T" + to_string(i + 1) + ".txt").c_str());
        }
        else
        {
            next = fileRead(i + 1, tsv);
            string str;
            getline(tfile, str);
            stringstream s(str);
            getline(s, str, ' ');
            int8_t save = stoi(str);
        }
        tfile.close();

        long saveEnd = current.size() * 15;
        while (save < saveEnd)
        {
            vector<thread> threads = {};
            for (int8_t i = 0; i < numThreads; i++)
            {
                threads.emplace_back(thread(threadMult, ref(threadVectors[i]), i, save, ref(tsv), ref(prior), ref(current)));
            }
            // Do matrix multiplication in threads
            for (int8_t i = 0; i < numThreads; i++)
            {
                threads[i].join();
                vector<SO6>::iterator itr = threadVectors[i].begin();
                vector<SO6>::iterator end = threadVectors[i].end();
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
    chrono::duration<double> timeelapsed = chrono::high_resolution_clock::now() - tbefore;
    std::cout << "\nTotal time elapsed: " << chrono::duration_cast<chrono::milliseconds>(timeelapsed).count() << "ms\n";
    return 0;
}
