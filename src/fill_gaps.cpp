
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

void run(NeighborFunction& NF, vector<string>& seqs, vector<string>& baits, vector<vector<bool>> & cover_marks, LL max_gap, LL bait_length, bool verbose){
    for(LL seq_id = 0; seq_id < seqs.size(); seq_id++){
        vector<bool>& v = cover_marks[seq_id];
        LL idx = 0;
        while(idx < v.size()){
            LL run = run_length(v,idx);
            if(verbose && run > max_gap && v[idx] == 0) cout << "Filling a gap of length " << run << endl;
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
        }
    }
}

int main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Fills gaps in a given bait cover.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("G,max-gap", "Maximum allowable gap length.", cxxopts::value<LL>()->default_value("0"))
      ("d,hamming-distance", "Number of allowed mismatches in the baits.", cxxopts::value<LL>()->default_value("40"))
      ("s,sequences", "Path to the fasta file of the input sequences.", cxxopts::value<string>()->default_value(""))
      ("b,baits", "Path to the fasta file of the baits.", cxxopts::value<string>()->default_value(""))
      ("c,cover-marks", "Path to the file of the cover marks created by the greedy algorithm.", cxxopts::value<string>()->default_value(""))
      ("o,out", "Output filename (fasta)", cxxopts::value<string>()->default_value(""))
      ("fm-index-out", "The algorithm is based on FM-index, which we build at the start. Building the index can take a lot of time and memory. Use this option to save the FM-index to disk so that you can later run the algorithm with different parameters re-using the same FM-index. (optional).", cxxopts::value<string>()->default_value(""))
      ("f,fm-index", "Path to a previously saved FM-index on disk (--fm-index-out). This option loads the FM index from disk instead of building it again.", cxxopts::value<string>()->default_value(""))
      ("t,n-threads", "Maximum number of parallel threads. The program is not very well optimized for parallel processing, so don't expect much of a speedup here.", cxxopts::value<LL>()->default_value("1"))
      ("g,seed-len", "Seed and extend g-mer seed length", cxxopts::value<LL>()->default_value("20"))
      ("v,verbose", "Print debug output", cxxopts::value<bool>()->default_value("false"))
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
    string fm_index_file = cli_params["fm-index"].as<string>();
    string fmi_out = cli_params["fm-index-out"].as<string>();
    LL n_threads = cli_params["n-threads"].as<LL>();
    omp_set_num_threads(n_threads);
    LL max_gap = cli_params["max-gap"].as<LL>();
    assert(max_gap > 0);
    string outfile = cli_params["out"].as<string>();
    bool verbose = cli_params["verbose"].as<bool>();
    check_writable(outfile);

    cerr << "Loading sequences" << endl;
    vector<string> seqs = read_sequences(cli_params["sequences"].as<string>());
    cerr << "Loading baits" << endl;
    vector<string> baits = read_sequences(cli_params["baits"].as<string>());
    cerr << "Loading cover marks" << endl;
    vector<vector<bool> > cover_marks = read_cover_marks(cli_params["cover-marks"].as<string>());
    
    assert(baits.size() > 0);
    for(LL i = 1; i < baits.size(); i++) assert(baits[i].size() == baits[i-1].size());
    
    LL bait_length = baits[0].size();

    FM_index fmi;
    if(fm_index_file == ""){
        cerr << "Constructing FM index" << endl;
        fmi.construct(seqs);
        cerr << "FM index construction done" << endl;
        if(fmi_out != ""){
            cerr << "Saving the FM-index to " << fmi_out << endl;
            throwing_ofstream fmi_out_stream(fmi_out, ios_base::binary);
            fmi.serialize(fmi_out_stream.stream);
        }
    } else{
        throwing_ifstream in(fm_index_file, ios_base::binary); // Throws on error
        cerr << "Loading FM-index" << endl;
        fmi.load(in.stream);
        cerr << "FM index loaded" << endl;
    }

    FM_NeighborCandidateFunction NCF;
    NCF.init(&fmi, g);
    NeighborFunction NF;
    NF.init(&NCF, &seqs, d, bait_length);

    cerr << "Running" << endl;   
    run(NF, seqs, baits, cover_marks, max_gap, bait_length, verbose);

    cerr << "Writing result to " << outfile << endl;
    throwing_ofstream out(outfile);
    for(LL i = 0; i < baits.size(); i++){
        out << ">" << i << "\n" << baits[i] << "\n";
    }

}    