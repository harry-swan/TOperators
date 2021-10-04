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
#include <mutex>
#include <utility>
#include <map>
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

void generatePermutations(vector<vector<vector<char>>> p, vector<vector<vector<char>>> n,
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

int main()
{
    auto startTime = chrono::high_resolution_clock::now();
    ifstream p("patterns.txt");
    char pattern[255];
    list<pair<vector<vector<vector<char>>>, string>> patterns;
    while (p.getline(pattern, 255))
    {
        pair<vector<vector<vector<char>>>, string> pr;
        pr.first = readPattern(pattern);
        string s(pattern);
        pr.second = s;
        patterns.emplace_back(pr);
    }
    p.close();
    list<pair<vector<vector<vector<char>>>, string>>::iterator itr = patterns.begin();
    while (itr != patterns.end())
    {
        string patternName;
        patternName = "data/" + itr->second.substr(itr->second.find("["));
        fstream outputPatterns(patternName, fstream::out | fstream::app);
        fstream inputPatterns("../data/" + itr->second, fstream::in);
        outputPatterns << inputPatterns.rdbuf();
        inputPatterns.close();
        list<vector<vector<vector<char>>>> perms = {};
        generatePermutations(itr->first, {}, perms);
        list<pair<vector<vector<vector<char>>>, string>>::iterator itr2 = next(itr);
        while (itr2 != patterns.end())
        {
            if (equals(perms, itr2->first))
            {
                inputPatterns = fstream("../data/" + itr2->second, fstream::in);
                outputPatterns << inputPatterns.rdbuf();
                itr2--;
                patterns.erase(next(itr2));
                inputPatterns.close();
            }
            itr2++;
        }
        outputPatterns.close();
        itr++;
    }
    auto endTime = chrono::high_resolution_clock::now();
    auto ret = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    cout << ">>>Combined row permutation equivalent patterns in " << ret << "ms\n";
}