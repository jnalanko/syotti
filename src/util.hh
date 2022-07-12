#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <set>
#include <map>
#include <sstream>
#include <fstream>
#include "input_reading.hh"
#include <unordered_set>
#include <cmath>
#include "stdlib_printing.hh"
#include "progress.hh"

typedef int64_t LL;

string get_rc(string S);
string get_canonical(const string& S);
vector<string> split(string s, char delimiter);

// Get smallest g-mer of S
string get_minimizer(const string& S, LL g);

// 'N' is defined to not match to anything
LL hamming_distance_not_matching_N(const string& A, const string& B);

vector<string> read_sequences(string filename, bool append_reverse_complements = false);

template<typename T>
void write_to_file(vector<T>& v, string filename){
    throwing_ofstream out(filename);
    for(auto x : v) out << x << "\n";
}

vector<string> get_all_distinct_kmers_including_rc(const vector<string>& v, LL k);
void check_writable(string filename);
