#include <iostream>
#include <vector>
#include <thread>
#include <future>
#include <fstream>
#include <chrono>
#include <list>
#include <string>
#include <sstream>
#include <memory>
#include <functional>
#include <utility>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "../Z2.hpp"
#include "../SO6.hpp"
#include "../T_Hist.hpp"

using namespace std;

const unsigned char numThreads = 4;

struct Pattern
{
    vector<vector<vector<char>>> p;
    string s;
    bool e;
    shared_ptr<struct Pattern> next;
};

shared_ptr<struct Pattern> head;

void insert(vector<vector<vector<char>>> p, string s)
{
    shared_ptr<struct Pattern> pattern(new struct Pattern);
    pattern->p = p;
    pattern->s = s;
    pattern->next = head;
    pattern->e = true;
    head = pattern;
}

void
generatePermutations(vector<vector<vector<char>>> p, vector<vector<vector<char>>> n,
                     list<vector<vector<vector<char>>>> &perms)
{
    if (!p.size())
    {
        perms.emplace_back(n);
        return;
    }
    for (unsigned int i = 0; i < p.size(); i++)
    {
        vector<vector<vector<char>>> tp = p;
        vector<vector<vector<char>>> tn = n;
        tn.emplace_back(p[i]);
        tp.erase(tp.begin() + i);
        generatePermutations(tp, tn, perms);
    }
    return;
}

bool equals(list<vector<vector<vector<char>>>> perms, vector<vector<vector<char>>> p)
{
    list<vector<vector<vector<char>>>>::iterator itr = perms.begin();
    while (itr != perms.end())
    {
        if (*itr == p)
            return true;
        itr++;
    }
    return false;
}

vector<vector<vector<char>>> readPattern(string pattern)
{
    vector<vector<vector<char>>> transposeMatrix = {};
    stringstream pat(pattern.substr(pattern.find("[")));
    for (int i = 0; i < 6; i++)
    {
        transposeMatrix.emplace_back(vector<vector<char>>());
        for (int j = 0; j < 6; j++)
        {
            transposeMatrix[i].emplace_back(vector<char>());
        }
    }
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                pat.get();
                char c;
                pat.get(c);
                //Subtract 30 because 0 is 30 in ASCII
                transposeMatrix[i][j].emplace_back(c - 30);
            }
        }
        pat.get();
    }
    return transposeMatrix;
}

void rowWrite(list<vector<vector<vector<char>>>> perms, string patternName, shared_ptr<struct Pattern> start, vector<unsigned int> &counters, unsigned int idx)
{
    shared_ptr<struct Pattern> pattern = start;
    fstream inputPatterns;
    fstream outputPatterns(patternName, fstream::out | fstream::app);
    while(pattern)
    {
        if(pattern->e && equals(perms, pattern->p))
        {
            pattern->e = false;
            counters[idx]++;
            inputPatterns = fstream("../data/" + pattern->s, fstream::in);
            outputPatterns << inputPatterns.rdbuf();
            inputPatterns.close();
        }
        pattern = pattern->next;
    }
}

int main()
{
    auto startTime = chrono::high_resolution_clock::now();
    ifstream p("patterns.txt");
    char pattern[255];
    unsigned int e_count = 0;
    while (p.getline(pattern, 255))
    {
        string s(pattern);
        insert(readPattern(pattern), s);
        e_count++;
    }
    p.close();
    list<vector<vector<vector<char>>>> threadPatterns;
    atomic<int> activeThreads(0);
    while (e_count > 0)
    {
        shared_ptr<struct Pattern> newPattern;
        newPattern = head;
        vector<thread> threads = {};
        vector<unsigned int> counters(numThreads);
        unsigned int counter_pos = 0;
        while(newPattern && activeThreads < numThreads)
        {
            if(newPattern->e)
            {
                list<vector<vector<vector<char>>>>::iterator itr = threadPatterns.begin();
                list<vector<vector<vector<char>>>> perms = {};
                generatePermutations(newPattern->p, {}, perms);
                while(itr != threadPatterns.end() && !equals(perms, *itr))
                {
                    itr++;
                }
                if(itr == threadPatterns.end())
                {
                    string patternName = "data/" + newPattern->s.substr(newPattern->s.find("["));
                    threadPatterns.emplace_back(newPattern->p);
                    threads.emplace_back(thread(rowWrite, perms, patternName, newPattern, ref(counters), counter_pos++));
                    activeThreads++;
                }
            }
            newPattern = newPattern->next;
        }
        unsigned short usedThreads = activeThreads;
        counter_pos = 0;
        while(activeThreads > 0)
        {
            threads[usedThreads - activeThreads].join();
            activeThreads--;
            e_count -= counters[counter_pos++];
        }
        threadPatterns.clear();
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto ret = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    std::cout << ">>>Combined row permutation equivalent patterns in " << ret << "ms\n";
}