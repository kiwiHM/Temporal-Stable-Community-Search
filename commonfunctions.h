#ifndef COMMONFUNCTIONS
#define COMMONFUNCTIONS

// #include <bits/stdc++.h>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <queue>
#include <list>
#include <assert.h>

using namespace std;

using pii = pair <int, int>;
using piii = pair <int, pii>;
using piiii = pair <pii, pii>;

void putProcess(double procedure, double time_used);

std::stringstream timeFormatting(double microSeconds);

unsigned long long currentTime();

bool checksame_vector(vector <int> vec1, vector <int> vec2);

// template <typename tn> void read(tn &a);
void read(int &a);

struct pair_hash{
    inline size_t operator()(const pii & v) const {
        return v.first * 32 + v.second;
    }
};

void print_vec(vector <int> vec, string st = "", bool flag = true);

#endif