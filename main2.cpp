/**
 * T Operator Product Generation Main File
 * @file main.cpp
 * @author Swan Klein
 * @author Connor Mooney
 * @author Michael Jarret
 * @author Andrew Glaudell
 * @author Jacob Weston
 * @author Mingzhen Tian
 * @version 8/15/21
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <thread>
#include <future>
#include <fstream>
#include <chrono>
#include <cmath>
#include <set>
#include <list>
#include <string>
#include <sstream>
#include <functional>
#include <mutex>
#include <shared_mutex>
#include <utility>
#include <map>
#include <atomic>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "Z2.hpp"
#include "SO6.hpp"
#include "T_Hist.hpp"

using namespace std;

const unsigned char numThreads = 4;
unsigned long long operationsPerThread;
unsigned char rem;

const unsigned char tCount = 4;

// For this and above, brute force results into a vector
// Has no effect if setless > tCount
const unsigned char setless = 9;

//Turn this on if you want to read in saved data
//Make sure you have the txt files for what genFrom is set to if true
const bool tIO = false;

//If tIO true, choose which tCount to begin generating from:
const unsigned char genFrom = tCount;

//Saves every saveInterval iterations
//This also determines parallel block sizes
unsigned long long saveInterval = 50000 * numThreads;

//A Map of patterns to their lock and file, plus a lock for appending to this map
map<string, pair<shared_ptr<mutex>, shared_ptr<vector<T_Hist>>>> patternMap;
shared_timed_mutex mapLock;
atomic<unsigned long> counter;
const unsigned long counterWrite = 10000000;

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

// File reading from a partially completed T list
// Defunct now that we are printing patterns
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
    //T_Hist m;
    while (tfile.get(hist))
    {
        //Convert hex character to integer
        tmp.push_back((hist >= 'a') ? (hist - 'a' + 10) : (hist - '0'));
        if (++i % tc == 0)
        {
            tset.insert(T_Hist(tmp));
            tmp.clear();
            //tfile.get(hist);
        }
    }
    return tset;
}

// This is used to write results to the files named after patterns
void writeResults(unsigned char i, unsigned char threadNum, list<T_Hist> &append)
{
    unsigned long long start = threadNum * (append.size() / numThreads) + min(threadNum, (unsigned char)(append.size() % numThreads));
    unsigned long long end = start + (append.size() / numThreads) + (threadNum < (unsigned char)(append.size() % numThreads));
    list<T_Hist>::iterator itr = append.begin();
    list<T_Hist>::iterator itrend = append.begin();
    advance(itr, start);
    advance(itrend, end);
    while (itr != itrend)
    {
        stringstream patternName;
        patternName << "data/" << itr->reconstruct() << ".txt";
        string name = patternName.str();
        const shared_lock<shared_timed_mutex> *slck = new shared_lock<shared_timed_mutex>(mapLock);
        map<string, pair<shared_ptr<mutex>, shared_ptr<vector<T_Hist>>>>::const_iterator mapItr = patternMap.find(name);
        bool found = (mapItr != patternMap.end());
        pair<shared_ptr<mutex>, shared_ptr<vector<T_Hist>>> patternPair = mapItr->second;
        delete slck;
        if (!found)
        {
            const lock_guard<shared_timed_mutex> lck(mapLock);
            if (patternMap.find(name) == patternMap.end())
            {
                shared_ptr<mutex> pairLock(new mutex());
                shared_ptr<vector<T_Hist>> pairVector(new vector<T_Hist>());
                patternPair = patternMap.emplace(name, pair<shared_ptr<mutex>, shared_ptr<vector<T_Hist>>>(pairLock, pairVector)).first->second;
            }
            else
            {
                mapItr = patternMap.find(name);
                patternPair = mapItr->second;
            }
        }
        shared_ptr<mutex> appendLock = patternPair.first;
        shared_ptr<vector<T_Hist>> tVector = patternPair.second;
        T_Hist hist = *itr++;
        const lock_guard<mutex> *lck = new lock_guard<mutex>(*appendLock);
        tVector->emplace_back(hist);
        delete lck;
        counter++;
    }
}

// Outputs the results in patternMap to files
void outputResults(unsigned char i, unsigned char threadNum)
{
    unsigned long long start = threadNum * (patternMap.size() / numThreads) + min(threadNum, (unsigned char)(patternMap.size() % numThreads));
    unsigned long long end = start + (patternMap.size() / numThreads) + (threadNum < (unsigned char)(patternMap.size() % numThreads));
    map<string, pair<shared_ptr<mutex>, shared_ptr<vector<T_Hist>>>>::iterator itr = patternMap.begin();
    map<string, pair<shared_ptr<mutex>, shared_ptr<vector<T_Hist>>>>::iterator itrend = patternMap.begin();
    advance(itr, start);
    advance(itrend, end);
    auto startTime = chrono::high_resolution_clock::now();
    while (itr != itrend)
    {
        string name = itr->first;
        vector<T_Hist>::iterator vectItr = itr->second.second->begin();
        fstream outputFile(name, fstream::out | fstream::app);
        vector<T_Hist>::iterator vectEnd = itr->second.second->end();
        while (vectItr != vectEnd)
        {
            outputFile << *vectItr++ << ' ';
        }
        itr++->second.second->clear();
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto ret = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    cout << ">>>Wrote " << (end - start) << " residue patterns in " << ret << "ms\n";
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

void threadMultSetless(vector<T_Hist> &threadVector, const unsigned char threadNum, const unsigned long long threadContinue,
                       vector<T_Hist> &currentvec)
{
    // Setting up thread iteration, thread multiplies everything in range (titr, citr) to (tend, cend)
    // by all 15 base T Operators
    unsigned long long start = threadContinue + threadNum * operationsPerThread + min(threadNum, rem);
    if (start >= currentvec.size() * 15)
        return;
    unsigned long long end = min(start + operationsPerThread + (threadNum < rem), (unsigned long long)currentvec.size() * 15);
    vector<T_Hist>::iterator citr = currentvec.begin();
    vector<T_Hist>::iterator cend = currentvec.begin();
    // Add 1 because we don't care about multiplication by identity tsv[0]
    unsigned char t = 1 + start / currentvec.size();
    advance(citr, start % currentvec.size());
    unsigned char tend = 1 + end / currentvec.size();
    advance(cend, end % currentvec.size());

    T_Hist m, curr;
    while (t < 16)
    {
        while (citr != currentvec.end())
        {
            if (t == tend && citr == cend)
                break;
            curr = *citr;
            vector<unsigned char> vec = {t};
            m = T_Hist(vec) * curr;
            threadVector.push_back(m);
            citr++;
        }
        if (t == tend)
            break;
        t++;
        citr = currentvec.begin();
    }
}

int main()
{
    counter = 0;
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
    unsigned long long nSize = 0;
    list<T_Hist> append({});
    ifstream tfile;
    unsigned char start = 0;

    // Initialize the head of the SO6 tree
    T_Hist::initHead();
    // Generate the lookup table
    std::cout << "\nBeginning Table Generation with Depth " << (+tCount + 1) / 2 << "\n";

    T_Hist::tableInsert(T_Hist::head, NULL, (unsigned char)((+tCount + 1) / 2));

    if (tIO && genFrom > 2)
    {
        prior = fileRead(genFrom - 2);
        current = fileRead(genFrom - 1);
        start = genFrom - 1;
    }

    // /* Here's a thing I'm toying with
    // *  >>>>>>>>>>>>>>>>>>>>>>>>
    // *  >>>>>>>>>>>>>>>>>>>>>>>>
    // *  >>>>>>>>>>>>>>>>>>>>>>>>
    // */

    // std::cout << "Starting Test 1:\n";
    // bool test = true;
    // for(int tttt = 0; tttt < 100000; tttt++) {
    //     std::vector<unsigned char> vec;
    //     for(int m = 0; m < 5; m++) {
    //         unsigned char res = 1 + (rand() %15);
    //         vec.push_back(res);
    //     }
    //     T_Hist tmp = T_Hist(vec);
    //     SO6 first = tmp.reconstruct();
    //     SO6 second = first;
    //     second.reduced_rep();
    //     std::vector<int8_t> perm = SO6::lexicographic_permutation(first);

    //     for(int row = 0; row<6; row++) {
    //         for(int col = 0; col <6; col ++) {
    //             int c = perm[col];
    //             if(c < 0) {test = -first[-c-1][row]==second[col][row];}
    //             else {test = first[c-1][row]==second[col][row];}
    //             if(!test) break;
    //         }
    //         if(!test) break;
    //     }
    //     if(!test) break;

    //     for(int row = 0; row<6; row++) {
    //             for(int col = 0; col <6; col ++) {
    //                 int c = perm[col];
    //                 bool sign = c <0;
    //                 test = !(SO6::lexLess(first[abs(c)-1],second[col],sign,0) || SO6::lexLess(second[col],first[abs(c)-1],0,sign)) ;
    //                 if(!test) break;
    //             }
    //             if(!test) break;
    //     }
    //     if(!test) break;
    // }

    // if(test) std::cout << "Successful.\n";
    // else {std::cout << "Failed.\n";}

    // Get every T operator

    long long timer[tCount];

    unsigned char i;
    for (i = start; i < tCount && i < (setless - 1); i++)
    {
        std::cout << "\nBeginning T-Count " << (i + 1) << "\n";

        auto start = chrono::high_resolution_clock::now();

        // Main loop here
        unsigned long long save = 0;
        tfile.open(("data/T" + to_string(i + 1) + ".sav").c_str());
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
            save = stoi(str);
        }
        tfile.close();

        unsigned long long saveEnd = current.size() * 15;
        while (save < saveEnd)
        {
            vector<thread> threads = {};
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads.emplace_back(thread(threadMult, ref(threadVectors[j]), j, save, ref(prior), ref(current)));
            }
            // Do matrix multiplication in threads
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads[j].join();
                vector<T_Hist>::iterator itr = threadVectors[j].begin();
                vector<T_Hist>::iterator end = threadVectors[j].end();
                while (itr != end)
                {
                    T_Hist &tmp = *itr;
                    T_Hist::curr_history = &tmp;
                    SO6 tmp2 = T_Hist::curr_history->reconstruct();
                    tmp2.reduced_rep();
                    T_Hist::curr = &tmp2;
                    // next.insert(*itr);
                    if (next.insert(*itr).second)
                    {
                        // This is only called if the insert into next succeeded
                        append.push_back(*itr);
                    }
                    itr++;
                }
                threadVectors[j].clear();
            }
            save += saveInterval;
            threads.clear();
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads.emplace_back(thread(writeResults, i, j, ref(append)));
            }
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads[j].join();
            }
            append.clear();
            threads.clear();
            if (counter > counterWrite)
            {
                for (unsigned char j = 0; j < numThreads; j++)
                {
                    threads.emplace_back(thread(outputResults, i, j));
                }
                for (unsigned char j = 0; j < numThreads; j++)
                {
                    threads[j].join();
                }
                counter = 0;
            }
        }

        // End main loop
        auto end = chrono::high_resolution_clock::now();
        // Begin reporting
        auto ret = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        nSize = next.size();
        std::cout << ">>>Found " << nSize << " new matrices in " << ret << "ms\n";
        timer[i] = ret;
        prior.clear();
        prior.swap(current); // T++
        current.swap(next);  // T++
    }

    vector<T_Hist> currentvec(current.size());
    vector<T_Hist> nextvec();
    if (i < tCount)
    {
        copy(current.begin(), current.end(), currentvec.begin());
        saveInterval = nSize;
        operationsPerThread = saveInterval / numThreads;
        rem = saveInterval % numThreads;
    }
    current.clear();
    prior.clear();

    while (i < tCount)
    {
        unsigned int sCount = pow(15, 2 + i - setless);
        std::cout << "\nBeginning T-Count " << (i + 1) << "\n";
        auto start = chrono::high_resolution_clock::now();

        for (unsigned int j = 0; j < sCount; j++)
        {
            if (j % 15 == 0 && i > tCount)
            {
                ifstream tfile;
                tfile.open(("data/T" + to_string(i + 1) + 's' + to_string(1 + (j / 15)) + ".txt").c_str());
                if (!tfile)
                {
                    cout << "File does not exist.\n";
                    exit(1);
                }
                currentvec.clear();
                char hist;
                unsigned long long k = 0;
                vector<unsigned char> tmp;
                //T_Hist m;
                while (tfile.get(hist))
                {
                    //Convert hex character to integer
                    tmp.push_back((hist >= 'a') ? (hist - 'a' + 10) : (hist - '0'));
                    if (++k % i == 0)
                    {
                        currentvec.push_back(T_Hist(tmp));
                        tmp.clear();
                    }
                }
            }
            vector<thread> threads = {};
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads.emplace_back(thread(threadMultSetless, ref(threadVectors[j]), j, 0, ref(currentvec)));
            }
            // Do matrix multiplication in threads
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads[j].join();
                vector<T_Hist>::iterator itr = threadVectors[j].begin();
                vector<T_Hist>::iterator end = threadVectors[j].end();
                while (itr != end)
                {
                    T_Hist &tmp = *itr;
                    T_Hist::curr_history = &tmp;
                    append.push_back(*itr);
                    itr++;
                }
                threadVectors[j].clear();
            }
            threads.clear();
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads.emplace_back(thread(writeResults, i, j, ref(append)));
            }
            for (unsigned char j = 0; j < numThreads; j++)
            {
                threads[j].join();
            }
            append.clear();
            threads.clear();
            if (counter > counterWrite)
            {
                for (unsigned char j = 0; j < numThreads; j++)
                {
                    threads.emplace_back(thread(outputResults, i, j));
                }
                for (unsigned char j = 0; j < numThreads; j++)
                {
                    threads[j].join();
                }
                counter = 0;
            }
        }

        // End main loop
        auto end = chrono::high_resolution_clock::now();
        // Begin reporting
        auto ret = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        nSize *= 15;
        std::cout << ">>>Found " << nSize << " new matrices in " << ret << "ms\n";
        timer[i] = ret;
        i++;
    }

    vector<thread> threads = {};
    for (unsigned char j = 0; j < numThreads; j++)
    {
        threads.emplace_back(thread(outputResults, i, j));
    }
    for (unsigned char j = 0; j < numThreads; j++)
    {
        threads[j].join();
    }

    // Free all table memory
    T_Hist::tableDelete(T_Hist::head, NULL);
    patternMap.clear();
    chrono::duration<double> timeelapsed = chrono::high_resolution_clock::now() - tbefore;
    std::cout << "\nTotal time elapsed: " << chrono::duration_cast<chrono::milliseconds>(timeelapsed).count() << "ms\n\n\n";
    std::cout << "{";
    for (int i = start; i < tCount - 1; i++)
    {
        std::cout << timer[i] << ",";
    }
    std::cout << timer[tCount - 1] << "}\n";
    return 0;
}