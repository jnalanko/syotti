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
#include "neighborhood_structure.hh"
#include "neighborhood_functions.hh"
#include "progress.hh"
#include "util.hh"
#include "greedy.hh"
#include "cxxopts.hpp"
#include <omp.h>

using namespace std;

bool is_branching(const string& bait, FM_index& fmi, LL dbg_order){
    LL bait_length = bait.size();
    FM_index::Interval I = fmi.search(bait.substr(bait_length-dbg_order, dbg_order));
    LL freq = I.size();
    for(LL i = bait_length-dbg_order-1; i >= 0; i--){
        I = fmi.left_extend(I, bait[i]);
        if(I.size() != freq) return true;
    }
    return false;
}

//    Backward search all the way from end
//        If the pattern of length k has a different frequency than pattern of length (k-1)
//            There is a branch-crossing at order (k-1)

// end-index is inclusive
vector<bool> detect_branch_crossing_orders(const string& bait, FM_index& fmi, LL end){
    vector<bool> ans(bait.size()+1); // ans[i] = 1 iff there is a branch detected at order i
    FM_index::Interval I = fmi.get_full_interval();
    LL prev_freq = -1;
    for(LL i = end; i >= 0; i--){
        I = fmi.left_extend(I, bait[i]);
        if(i < end && I.size() != prev_freq){
            ans[end - i] = true;
        }
        prev_freq = I.size();
    }
    return ans;
}

// Returns vector v such that v[i] = number of baits hitting a branch crossing at DBG order i
vector<LL> count_branch_crossings(vector<string>& sequences, vector<string>& baits, FM_index& fmi, LL n_threads){
    
    LL bait_length = baits[0].size();
    Progress_printer pp(baits.size(), 100);
    vector<LL> crossing_counts(bait_length+1);
    
    for(LL bait_id = 0; bait_id < baits.size(); bait_id++){
        string bait = baits[bait_id];
        string bait_rc = get_rc(bait);
        vector<bool> crossing_orders(bait.size() + 1);
        #pragma omp parallel for
        for(LL end = 1; end < bait.size(); end++){
            vector<bool> fw = detect_branch_crossing_orders(bait, fmi, end);
            vector<bool> bw = detect_branch_crossing_orders(bait_rc, fmi, end);
            assert(fw.size() == bw.size() && crossing_orders.size() == fw.size());
            #pragma omp critical
            for(LL i = 0; i < fw.size(); i++) crossing_orders[i] = crossing_orders[i] | fw[i] | bw[i];
        }
        for(LL order = 1; order < bait.size(); order++)
            crossing_counts[order] += crossing_orders[order];
        pp.job_done();
    }

    return crossing_counts;

}

// Returns a sorted range of integers including the endpoints such that the rest of the points are evenly
// distributed in between
vector<LL> linspace(LL start, LL end, LL n_points){
    assert(n_points >= 2);
    LL step = (end-start+1)/(n_points - 2 + 1);
    vector<LL> ans = {start};
    for(LL i = 0; i < n_points-2; i++)
        ans.push_back(ans.back() + step);
    ans.push_back(end);
    return ans;
}

// Iteration state of the function run. Passed to callbacks.
struct State{
    LL total_to_cover = 0;
    LL n_covered = 0;
    LL bait_id = 0;
    string bait;
    vector<vector<bool> > cover;
};

typedef std::function<void(const State&)> callback_t; // Function taking the bait_id of the latest bait and the current cover

void run(vector<string>& seqs, vector<string>& baits, FM_index& fmi, LL d, LL g, LL bait_length, bool verbose, vector<callback_t>& callbacks){
    State state;

    FM_NeighborCandidateFunction NCF;
    NCF.init(&fmi, g);
    NeighborFunction NF;
    NF.init(&NCF, &seqs, d, bait_length);
    for(string& S : seqs){
        vector<bool> v(S.size());
        state.total_to_cover += S.size();
        state.cover.push_back(v);
    }

    //vector<LL> callback_indices = linspace(0, baits.size()-1, n_callback_points);
    //std::reverse(callback_indices.begin(), callback_indices.end());
    Progress_printer pp(baits.size(), 100);
    
    state.n_covered = 0;
    for(state.bait_id = 0; state.bait_id < baits.size(); state.bait_id++){
        state.bait = baits[state.bait_id];
        vector<pair<LL,LL> > neighbors = NF.get_neighbors(state.bait);
        for(pair<LL,LL> P : NF.get_neighbors(get_rc(state.bait))) neighbors.push_back(P);
        for(pair<LL,LL> P : neighbors){
            LL doc_id, pos; std::tie(doc_id, pos) = P;
            if(verbose){
                cout << "--" << "\n" << state.bait << "\n" << seqs[doc_id].substr(pos, bait_length) << "\n";
            }
            assert(hamming_distance(state.bait, seqs[doc_id].substr(pos, bait_length)) <= d ||
                   hamming_distance(get_rc(state.bait), seqs[doc_id].substr(pos, bait_length)) <= d);
            for(LL i = pos; i < pos+bait_length; i++){
                if(state.cover[doc_id][i] == 0) state.n_covered++;
                state.cover[doc_id][i] = 1;
            }            
        }
        pp.job_done(to_string((double)state.n_covered/state.total_to_cover));
        for(auto& callback : callbacks) callback(state);
        /*if(callback_indices.size() > 0 && callback_indices.back() == state.bait_id){
            cerr << "Running callbacks" << endl;
            for(auto& callback : callbacks) callback(state);
            callback_indices.pop_back();
        }*/
    }
}

template<typename T>
vector<T> read_things_from_file(string filename){
    T x;
    throwing_ifstream in(filename);
    vector<T> ans;
    while(in.stream >> x) ans.push_back(x);
    return ans;
}

template<typename T>
LL run_length(const vector<T>& v, LL from){
    LL idx = from;
    while(idx < v.size() && v[idx] == v[from]) idx++;
    return idx-from;
}

int main(int argc, char** argv){

    cxxopts::Options options("Examine", "Computes various statistics");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("v,verbose", "Print more debug output", cxxopts::value<bool>()->default_value("false"))
      ("o,out-prefix", "The prefix for the output files", cxxopts::value<string>()->default_value(""))
      ("f,fm-index", "Path to the FM index file of the sequences (optional)", cxxopts::value<string>()->default_value(""))
      ("s,sequences", "Path to the fasta file of the input sequences", cxxopts::value<string>()->default_value(""))
      ("t,n-threads", "Number of parallel threads", cxxopts::value<LL>()->default_value("1"))
      ("b,baits", "Path to the fasta file of the baits", cxxopts::value<string>()->default_value(""))
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

    LL g = cli_params["g"].as<LL>();
    LL d = cli_params["d"].as<LL>();
    LL bait_length = cli_params["bait-len"].as<LL>();
    string out_prefix = cli_params["out-prefix"].as<string>();
    string baitfile = cli_params["baits"].as<string>();
    string fm_index_file = cli_params["fm-index"].as<string>();
    string sequence_file = cli_params["sequences"].as<string>();
    bool verbose = cli_params["verbose"].as<bool>();
    //bool compute_gaps = cli_params["gaps"].as<bool>();
    LL n_threads = cli_params["n-threads"].as<LL>();
    //LL n_points = cli_params["n-eval-points"].as<LL>();
    omp_set_num_threads(n_threads);

    if(fm_index_file == ""){
        cerr << "Error: FM index file not given" << endl; exit(1);
    }

    if(sequence_file == ""){
        cerr << "Error: Sequence file not given" << endl; exit(1);
    }

    if(out_prefix == ""){
        cerr << "Error: Prefix not given" << endl; exit(1);
    }

    if(baitfile == ""){
        cerr << "Error: Bait file not given" << endl; exit(1);
    }

    cerr << "Loading sequences" << endl;
    vector<string> sequences = read_sequences(sequence_file, true); // Reverse complements included to match the FM index
    cerr << "Loading baits" << endl;
    vector<string> baits = read_sequences(baitfile, false); // Reverse complements not included

    cerr << "Loading FM index" << endl;
    FM_index fmi;
    throwing_ifstream fmi_in(fm_index_file);
    fmi.load(fmi_in.stream);

    throwing_ofstream final_gaps_out(out_prefix + "-gaps.txt");
    throwing_ofstream final_crossings_out(out_prefix + "-crossings.txt");
    throwing_ofstream final_cover_marks_out(out_prefix + "-cover_marks.txt");
    throwing_ofstream cover_curve_out(out_prefix + "-cover-curve.txt");

    callback_t gap_length_callback = [&](const State& state){
        if(state.bait_id == baits.size()-1){
            for(auto & v : state.cover){
                LL idx = 0;
                while(idx < v.size()){
                    if(v[idx] == 0) final_gaps_out.stream << run_length(v, idx) << " ";
                    idx += run_length(v,idx);
                }
            }
            final_gaps_out.stream << endl;
        }
    };

    callback_t write_cover_curve_callback = [&](const State& state){
        cover_curve_out.stream << (double)state.n_covered/state.total_to_cover << "\n";
    };

    callback_t write_cover_marks_callback = [&](const State& state){
        
        if(state.bait_id == baits.size()-1){
            cerr << "Writing cover marks" << endl;
            for(const vector<bool>& v : state.cover){
                for(bool b : v) final_cover_marks_out.stream << b;
                final_cover_marks_out.stream << "\n";
            }
        }
    };

    callback_t count_branch_crossings_callback = [&](const State& state){
        if(state.bait_id == baits.size()-1){
            cerr << "Counting branch crossings" << endl;
            vector<LL> crossing_counts = count_branch_crossings(sequences, baits, fmi, n_threads);
            for(LL order = 1; order < crossing_counts.size(); order++){
                final_crossings_out.stream << order << " " << (double)crossing_counts[order]/baits.size() << endl;
            }
        }
    };

    vector<callback_t> callbacks = {gap_length_callback, write_cover_marks_callback, write_cover_curve_callback, count_branch_crossings_callback};

    run(sequences, baits, fmi, d, g, bait_length, verbose, callbacks);


}