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

string get_rc(string S){
    std::reverse(S.begin(), S.end());
    for(char& c : S){
        if(c == 'A') c = 'T';
        else if(c == 'C') c = 'G';
        else if(c == 'G') c = 'C';
        else if(c == 'T') c = 'A';
        else{
            cerr << "Error: invalid character in sequence: " << c << endl;
            cerr << "(Should be one of upper case characters A,C,G,T)" << endl;
            exit(1);
        }
    }
    return S;
}

string get_canonical(const string& S){
    string T = get_rc(S);
    return S < T ? S : T;
}


vector<string> split(string s, char delimiter){
    stringstream test(s);
    string segment;
    vector<string> seglist;

    while(getline(test, segment, delimiter)){
        seglist.push_back(segment);
    }
    return seglist;
}

// Get smallest g-mer of S
string get_minimizer(const string& S, LL g){
    assert(S.size() >= g);
    string M = get_canonical(S.substr(0,g));
    for(LL i = 1; i < S.size()-g+1; i++){
        string x = get_canonical(S.substr(i,g));
        if(x < M) M = x;
        
    }
    return M;
}

LL hamming_distance(const string& A, const string& B){
    assert(A.size() == B.size());
    LL ans = 0;
    for(LL i = 0; i < A.size(); i++){
        if(A[i] != B[i]) ans++;
    }
    return ans;
}

vector<string> read_sequences(string filename, bool append_reverse_complements = false){
    vector<string> seqs;
    Sequence_Reader sr(filename, FASTA_MODE);
    while(!sr.done()){
        seqs.push_back(sr.get_next_query_stream().get_all());
    }
    if(append_reverse_complements){
        LL n = seqs.size();
        for(LL i = 0; i < n; i++)
            seqs.push_back(get_rc(seqs[i]));
    }
    return seqs;
}


template<typename T>
void write_to_file(vector<T>& v, string filename){
    throwing_ofstream out(filename);
    for(auto x : v) out << x << "\n";
}

vector<string> get_all_distinct_kmers_including_rc(const vector<string>& v, LL k){
    cerr << "Listing distinct k-mers (rc included)" << endl;
    Progress_printer pp(v.size(), 100);
    unordered_set<string> hashset;
    for(string S : v){
        for(LL i = 0; i < S.size()-k+1; i++){
            hashset.insert(S.substr(i,k));
            hashset.insert(get_rc(S.substr(i,k)));
        }
        pp.job_done(to_string(k) + "-mers");
    }
    vector<string> ans;
    for(string S : hashset) ans.push_back(S);
    return ans;
}

void check_writable(string filename){
    ofstream F(filename);
    if(!F.good()){
        cerr << "Error writing to file: " << filename << endl;
        exit(1);
    }
}
