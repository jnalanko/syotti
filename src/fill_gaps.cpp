
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
#include "neighborhood_functions.hh"
#include "progress.hh"
#include "util.hh"
#include "cxxopts.hpp"
#include <omp.h>
#include <time.h> 



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
    bool randomize;

    Greedy(){}

    // Sequences must be the same and in the same order as what the neighbor candiate function was built for
    // NCF pointer does not transfer ownership.
    void init(NeighborFunction* NF, vector<string>* seqs, LL bait_length, LL d, LL g, bool randomize){
        this->bait_length = bait_length;
        this->d = d;
        this->g = g;
        this->NF = NF;
        this->seqs = seqs;
        this->randomize = randomize;
    }

    struct Result{
        vector<string> baits;
        vector<double> cover_fractions;
        // Taking baits[0..i] covers cover_fractions[i] of the bases.

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
        vector<vector<bool> > covered; // Covered[i][j] = number of times pos j of sequence i is covered
        LL total_to_cover = 0;
        LL n_covered = 0;
        for(string S : *seqs){
            vector<bool> V(S.size(), 0);
            covered.push_back(V);
            total_to_cover += S.size();
        }

        // Run covering algorithm
        Progress_printer pp2(seqs->size(), 100);
        LL n_good_candidates_found = 0;
        LL n_seqs_processed = 0;
        for(LL sequence_id : create_permutation(seqs->size(), randomize)){
            string S = seqs->at(sequence_id);
            for(LL i = 0; i < S.size(); i++){
                if(!covered[sequence_id][i]){
                    LL start = min(i, (LL)S.size() - bait_length);
                    string x = S.substr(start,bait_length);
                    for(LL off = 0; off < bait_length; off++){
                        if(covered[sequence_id][start+off] == 0) n_covered++;
                        covered[sequence_id][start+off] = 1;
                    }
                    vector<pair<LL,LL> > neighbors = NF->get_neighbors(x);
                    for(auto P : NF->get_neighbors(get_rc(x))) neighbors.push_back(P);
                    for(pair<LL,LL> P : neighbors){
                        LL doc_id = P.first;
                        LL pos = P.second;
                        string y = seqs->at(doc_id).substr(pos, bait_length);
                        n_good_candidates_found++;
                        for(LL off = 0; off < bait_length; off++){
                            if(covered[doc_id][pos+off] == 0) n_covered++;
                            covered[doc_id][pos+off] = 1;
                        }
                    }
                    
                    result.baits.push_back(x);
                    result.cover_fractions.push_back((double)n_covered / total_to_cover);
                }
            }

            pp2.job_done(to_string(result.baits.size()));
            if(n_seqs_processed % 1000 == 0) cerr << (double)n_covered / total_to_cover << endl;
            n_seqs_processed++;
        }

        cout << "Good candidates: " << n_good_candidates_found << endl;
        cout << "Bait set size " << result.baits.size() << endl;

        return result;
    }

};

vector<vector<bool> > read_cover_marks(string filename){
    throwing_ifstream in(filename);
    string line;
    vector<vector<bool> > v;
    while(getline(in.stream, line)){
        vector<bool> row;
        for(char c : line){
            assert(c == '0' or c == '1');
            row.push_back(c == '1' ? 1 : 0);
        }
        v.push_back(row);
    }
    return v;
}

template<typename T>
LL run_length(const vector<T>& v, LL from){
    LL idx = from;
    while(idx < v.size() && v[idx] == v[from]) idx++;
    return idx-from;
}

void run(NeighborFunction& NF, vector<string>& seqs, vector<string>& baits, vector<vector<bool>> & cover_marks, LL max_gap, LL bait_length){
    for(LL seq_id = 0; seq_id < seqs.size(); seq_id++){
        vector<bool>& v = cover_marks[seq_id];
        LL idx = 0;
        while(idx < v.size()){
            LL run = run_length(v,idx);
            //cout << seq_id << " " << idx << " " << run << " " << v[idx] << endl;
            if(run > max_gap && v[idx] == 0) cout << run << endl;
            if(v[idx] == 1 || (v[idx] == 0 && run <= max_gap)) idx += run;
            else{ // Select bait
                idx += max_gap;
                LL start = min(idx, (LL)seqs[seq_id].size() - bait_length);
                string x = seqs[seq_id].substr(start,bait_length);
                baits.push_back(x);
                vector<pair<LL,LL> > neighbors = NF.get_neighbors(x);
                for(auto P : NF.get_neighbors(get_rc(x))) neighbors.push_back(P);
                for(pair<LL,LL> P : neighbors){
                    LL doc_id = P.first;
                    LL pos = P.second;
                    for(LL off = 0; off < bait_length; off++){
                        cover_marks[doc_id][pos+off] = 1;
                    }
                }
            }
            /*
                string bait = seqs[seq_id].substr(idx, bait_length);
                while(bait.size() < bait_length) bait += bait.back(); // Pad
                baits.push_back(bait);
                for(LL i = idx; i < min(idx + bait_length, (LL)v.size()); i++){
                    v[i] = 1;
                }
            }*/
        }
    }
}

int main(int argc, char** argv){

    cxxopts::Options options("Greedy","Computes a greedy bait cover. A precomputed FM-index can be given.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("f,fm-index", "Path to the FM index file of the sequences", cxxopts::value<string>()->default_value(""))
      ("s,sequences", "Path to a fasta FILE of the input sequences", cxxopts::value<string>()->default_value(""))
      ("b,baits", "Path to a fasta file of the input baits", cxxopts::value<string>()->default_value(""))
      ("c,cover-marks", "Path to a fasta file of the cover marks", cxxopts::value<string>()->default_value(""))
      ("o,out", "Filename for the full new bait set", cxxopts::value<string>()->default_value(""))
      ("G,max-gap", "Maximum allowable gap length", cxxopts::value<LL>()->default_value("0"))
      ("t,n-threads", "Number of parallel threads", cxxopts::value<LL>()->default_value("1"))
      ("g,seed-len", "Seed and extend g-mer seed length", cxxopts::value<LL>()->default_value("20"))
      ("L,bait-len", "Length of the baits", cxxopts::value<LL>()->default_value("120"))
      ("d,hamming-distance", "Number of allowed mismatches", cxxopts::value<LL>()->default_value("40"))
      ("h,help", "Print instructions.", cxxopts::value<bool>()->default_value("false"))
    ;

    cxxopts::ParseResult cli_params = options.parse(argc, argv);
    if (cli_params["help"].as<bool>() == true || original_argc == 1){
        cerr << options.help() << endl;
        return 1;
    }

    if(cli_params["out"].as<string>() == ""){
        cerr << "Error: output file not specified" << endl;
        return 1;
    }

    LL g = cli_params["g"].as<LL>();
    LL d = cli_params["d"].as<LL>();
    LL bait_length = cli_params["bait-len"].as<LL>();
    LL n_threads = cli_params["n-threads"].as<LL>();
    omp_set_num_threads(n_threads);
    LL max_gap = cli_params["max-gap"].as<LL>();
    assert(max_gap > 0);
    string outfile = cli_params["out"].as<string>();
    check_writable(outfile);

    cerr << "Loading sequences" << endl;
    vector<string> seqs = read_sequences(cli_params["sequences"].as<string>(), true); // Also appends reverse complements
    cerr << "Loading baits" << endl;
    vector<string> baits = read_sequences(cli_params["baits"].as<string>(), false); // Don't append reverse complements
    cerr << "Loading cover marks" << endl;
    vector<vector<bool> > cover_marks = read_cover_marks(cli_params["cover-marks"].as<string>());
    
    cerr << "Loading FM-index" << endl;
    FM_index fmi;
    throwing_ifstream fmi_in(cli_params["fm-index"].as<string>(), ios_base::binary); // Throws on error
    fmi.load(fmi_in.stream);
    cerr << "FM index loaded" << endl;

    FM_NeighborCandidateFunction NCF;
    NCF.init(&fmi, g);
    NeighborFunction NF;
    NF.init(&NCF, &seqs, d, bait_length);

    cerr << "Running" << endl;   
    run(NF, seqs, baits, cover_marks, max_gap, bait_length);

    cerr << "Writing result to " << outfile << endl;
    throwing_ofstream out(outfile);
    for(string& S : baits) out << ">\n" << S << "\n";

}    