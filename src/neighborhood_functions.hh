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
#include "FM_index.hh"

/* THEFUNCTIONS IN THIS FILE SHOULD NOT KNOW ABOUT REVERSE COMPLEMENTS */

typedef int64_t LL;


// Abstract interface for the neighbor function
class NeighborCandidateFunction{

public:

    // Returns a pair (document id, starting position of a baitmer match in the document)
    virtual vector<pair<LL,LL> > get_candidates(const string& baitmer) = 0;

    virtual ~NeighborCandidateFunction() {}

};

class FM_NeighborCandidateFunction : public NeighborCandidateFunction{
private:

    FM_index* FMI; // Non-owning pointer
    LL g;

    // No copying because of pointer member
    FM_NeighborCandidateFunction(const FM_NeighborCandidateFunction& other) = delete;
    FM_NeighborCandidateFunction& operator=(const FM_NeighborCandidateFunction& other) = delete;

public:

    FM_NeighborCandidateFunction() {}


    // The class does not assume ownership of the FM-index. It will not de-allocate it.
    void init(FM_index* fmi, LL g){
        this->g = g;
        this->FMI = fmi;
    }

    // Returns a pair (document id, starting position of a baitmer match in the document)
    // This does not return reverse complement matches.
    vector<pair<LL,LL> > get_candidates(const string& baitmer){
        vector<vector<pair<LL,LL> > > ans_vecs(baitmer.size()-g+1);

        #pragma omp parallel for
        for(LL i = 0; i < baitmer.size()-g+1; i++){
            string z = baitmer.substr(i,g);
            ans_vecs[i] = FMI->subpattern_search(z, baitmer.size(), i);
        }

        // Collect answer vectors
        vector<pair<LL,LL> >  ans;
        for(auto& v : ans_vecs) for(auto x : v) ans.push_back(x);

        // Remove duplicates
        std::sort(ans.begin(), ans.end());
        ans.erase(std::unique(ans.begin(), ans.end() ), ans.end());
        return ans;
    }

};

class HashNeighborCandidateFunction : public NeighborCandidateFunction{

private:

    // g-mer x -> pairs (document id, position)
    unordered_map<string,vector<pair<LL,LL >>> gmer_map;
    vector<string> seqs;
    LL g;

public:



    // Baitmers: forward and rc versions of all baitmers
    HashNeighborCandidateFunction(){}

    void init(const vector<string>& seqs, LL bait_length, LL g){
        this->g = g;
        this->seqs = seqs;
        
        gmer_map.clear();
        cerr << "Computing g-mer map" << endl;
        Progress_printer pp1(seqs.size(), 100);
        for(LL seq_id = 0; seq_id < seqs.size(); seq_id++){
            string x = seqs[seq_id];
            for(LL i = 0; i < x.size()-g+1; i++){
                string z = x.substr(i,g);
                gmer_map[z].push_back({seq_id, i});
            }
            pp1.job_done("g-mers");
        }
    }

    vector<pair<LL,LL>> get_candidates(const string& baitmer){
        set<pair<LL,LL>> ans_set;
        LL L = baitmer.size();

        for(LL i = 0; i < baitmer.size()-g+1; i++){
            string z = baitmer.substr(i,g);
            for(pair<LL,LL> P : gmer_map[z]){
                LL seq_id = P.first;
                LL p = P.second - i;
                if(p >= 0 && p + L <= seqs[seq_id].size())
                ans_set.insert({seq_id, p});
            }
        }
        vector<pair<LL,LL>> ans(ans_set.begin(), ans_set.end());
        return ans;
    }

};


class NeighborFunction{

private:

    NeighborCandidateFunction* NCF; // Non-owning pointer
    vector<string>* seqs; // Non-owning pointer
    LL bait_length;
    LL d;


    NeighborFunction(const NeighborFunction& other) = delete;
    NeighborFunction& operator=(const NeighborFunction& other) = delete; 

public:

    NeighborFunction() {}

    // g = seed length, d = hamming distance bound
    void init(NeighborCandidateFunction* NCF, vector<string>* seqs, LL d, LL bait_length){
        this->NCF = NCF;
        this->seqs = seqs;
        this->bait_length = bait_length;
        this->d = d;
    }

    vector<pair<LL,LL> > get_neighbors(const string& baitmer){
        vector<pair<LL,LL> > candidates = NCF->get_candidates(baitmer);
        vector<pair<LL,LL> > verified_matches;
        for(pair<LL,LL> P : candidates){
            LL doc_id = P.first;
            LL pos = P.second;
            string y = seqs->at(doc_id).substr(pos, bait_length);
            if(hamming_distance_not_matching_N(baitmer,y) <= d)
                verified_matches.push_back({doc_id, pos});
        }
        return verified_matches;
    }

};
