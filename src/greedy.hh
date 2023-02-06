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
#include "SeqIO.hh"
#include <unordered_set>
#include <cmath>
#include "stdlib_printing.hh"
#include "neighborhood_functions.hh"
#include "progress.hh"
#include "util.hh"
#include <time.h> 

/*
 Algorithm sketch
 For each input string S:
     pos = 0
     If pos is covered 
         increment pos
     else 
        take the baitmer starting at pos. Mark all positions in all strings that are covered by this baitmer.
        Slide forward bait length.
*/

class Greedy{

friend class Greedy_tester;

private:

    // No copying this class because we have a pointer member and I don't want to deal with the copy semantics
    Greedy(const Greedy& other) = delete;
    Greedy& operator=(const Greedy& other) = delete;
    

public:
    LL bait_length; // Bait length
    LL d; // Number of allowed mismatches
    LL g; // Minimizer length
    vector<string>* seqs; // Non-owning pointer
    NeighborFunction* NF; // Non-owning pointer
    double cutoff; // Stop after this fraction has been covered
    bool randomize;

    Greedy(){}

    // Sequences must be the same and in the same order as what the neighbor candiate function was built for
    // NCF pointer does not transfer ownership.
    void init(NeighborFunction* NF, vector<string>* seqs, LL bait_length, LL d, LL g, bool randomize, double cutoff){
        this->bait_length = bait_length;
        this->d = d;
        this->g = g;
        this->NF = NF;
        this->seqs = seqs;
        this->randomize = randomize;
        this->cutoff = cutoff;
    }

    struct Result{
        vector<string> baits; // The final baits

        vector<double> cover_fractions; // Taking baits[0..i] covers cover_fractions[i] of the bases.
        
        vector<vector<bool> > covered; // cover[i][j] = 1 iff position j of sequence i is covered in the final cover

        Result(){}
    };

    vector<LL> create_permutation(LL n_elements, bool shuffle){
        vector<LL> v;
        v.resize(n_elements);
        for(LL i = 0; i < n_elements; i++) v[i] = i;
        if(shuffle){
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));
        }
        return v;
    }

    Result run(){
        Result result;
        // Initialize covering algorithm
        LL total_to_cover = 0;
        LL n_covered = 0;
        for(string S : *seqs){
            vector<bool> V(S.size(), 0);
            result.covered.push_back(V);
            total_to_cover += S.size();
        }

        // Run covering algorithm
        Progress_printer pp2(seqs->size(), 100);
        LL n_good_candidates_found = 0;
        LL n_seqs_processed = 0;
        bool done = false;
        for(LL sequence_id : create_permutation(seqs->size(), randomize)){
            string S = seqs->at(sequence_id);
            for(LL i = 0; i < S.size(); i++){
                if(!result.covered[sequence_id][i]){
                    LL start = min(i, (LL)S.size() - bait_length);
                    string x = S.substr(start,bait_length);
                    for(LL off = 0; off < bait_length; off++){
                        if(result.covered[sequence_id][start+off] == 0) n_covered++;
                        result.covered[sequence_id][start+off] = 1;
                    }
                    vector<pair<LL,LL> > neighbors = NF->get_neighbors(x);
                    for(auto P : NF->get_neighbors(get_rc(x))) neighbors.push_back(P);
                    for(pair<LL,LL> P : neighbors){
                        LL doc_id = P.first;
                        LL pos = P.second;
                        n_good_candidates_found++;
                        for(LL off = 0; off < bait_length; off++){
                            if(result.covered[doc_id][pos+off] == 0) n_covered++;
                            result.covered[doc_id][pos+off] = 1;
                        }
                    }
                    
                    result.baits.push_back(x);
                    result.cover_fractions.push_back((double)n_covered / total_to_cover);
                    if((double)n_covered / total_to_cover > cutoff){
                        done = true;
                        break;
                    }
                }
            }

            pp2.job_done(to_string(result.baits.size()));
            if(n_seqs_processed % 1000 == 0) cerr << (double)n_covered / total_to_cover << endl;
            n_seqs_processed++;
            if(done) break;
        }

        cout << "Good candidates: " << n_good_candidates_found << endl;
        cout << "Bait set size " << result.baits.size() << endl;

        return result;
    }

};