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
    vector<vector<int32_t> > cover; // Number of times each base pair is covered
    vector<pair<LL,LL> > neighbors;
};

typedef std::function<void(const State&)> callback_t; // Function taking the bait_id of the latest bait and the current cover

void run(vector<string>& seqs, vector<string>& baits, FM_index& fmi, LL d, LL g, LL bait_length, bool verbose, vector<callback_t>& callbacks){
    State state;

    FM_NeighborCandidateFunction NCF;
    NCF.init(&fmi, g);
    NeighborFunction NF;
    NF.init(&NCF, &seqs, d, bait_length);
    for(string& S : seqs){
        vector<int32_t> v(S.size());
        state.total_to_cover += S.size();
        state.cover.push_back(v);
    }

    //vector<LL> callback_indices = linspace(0, baits.size()-1, n_callback_points);
    //std::reverse(callback_indices.begin(), callback_indices.end());
    Progress_printer pp(baits.size(), 100);
    
    state.n_covered = 0;
    for(state.bait_id = 0; state.bait_id < baits.size(); state.bait_id++){
        state.bait = baits[state.bait_id];
        state.neighbors = NF.get_neighbors(state.bait);
        for(pair<LL,LL> P : NF.get_neighbors(get_rc(state.bait))) state.neighbors.push_back(P);
        for(pair<LL,LL> P : state.neighbors){
            LL doc_id, pos; std::tie(doc_id, pos) = P;
            if(verbose){
                cout << "--" << "\n" << state.bait << "\n" << seqs[doc_id].substr(pos, bait_length) << "\n";
            }
            assert(hamming_distance_not_matching_N(state.bait, seqs[doc_id].substr(pos, bait_length)) <= d ||
                   hamming_distance_not_matching_N(get_rc(state.bait), seqs[doc_id].substr(pos, bait_length)) <= d);
            for(LL i = pos; i < pos+bait_length; i++){
                if(state.cover[doc_id][i] == 0) state.n_covered++;
                state.cover[doc_id][i]++;
            }
        }
        pp.job_done(to_string((double)state.n_covered/state.total_to_cover));
        for(auto& callback : callbacks) callback(state);
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

    cxxopts::Options options(argv[0], "Computes various statistics on a given bait set.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("d,hamming-distance", "Number of allowed mismatches in the baits.", cxxopts::value<LL>()->default_value("40"))
      ("s,sequences", "Path to the fasta file of the input sequences.", cxxopts::value<string>()->default_value(""))
      ("b,baits", "Path to the fasta file of the baits.", cxxopts::value<string>()->default_value(""))
      ("o,out-prefix", "Filename prefix for the output files.", cxxopts::value<string>()->default_value(""))
      ("fm-index-out", "The algorithm is based on FM-index, which we build at the start. Building the index can take a lot of time and memory. Use this option to save the FM-index to disk so that you can later run the algorithm with different parameters re-using the same FM-index. (optional).", cxxopts::value<string>()->default_value(""))
      ("f,fm-index", "Path to a previously saved FM-index on disk (--fm-index-out). This option loads the FM index from disk instead of building it again.", cxxopts::value<string>()->default_value(""))
      ("t,n-threads", "Maximum number of parallel threads. The program is not very well optimized for parallel processing, so don't expect much of a speedup here.", cxxopts::value<LL>()->default_value("1"))
      ("g,seed-len", "The length of the seeds in the FM-index seed-and-extend approximate string search subroutine. A lower value will find more matches, but will be slower.", cxxopts::value<LL>()->default_value("20"))
      ("v,verbose", "Print debug output", cxxopts::value<bool>()->default_value("false"))
      ("h,help", "Print instructions.", cxxopts::value<bool>()->default_value("false"))
    ;

    cxxopts::ParseResult cli_params = options.parse(argc, argv);
    if (cli_params["help"].as<bool>() == true || original_argc == 1){
        cerr << options.help() << endl;
        return 1;
    }

    LL g = cli_params["g"].as<LL>();
    LL d = cli_params["d"].as<LL>();
    string out_prefix = cli_params["out-prefix"].as<string>();
    string baitfile = cli_params["baits"].as<string>();
    string fm_index_file = cli_params["fm-index"].as<string>();
    string fmi_out = cli_params["fm-index-out"].as<string>();
    string sequence_file = cli_params["sequences"].as<string>();
    bool verbose = cli_params["verbose"].as<bool>();
    LL n_threads = cli_params["n-threads"].as<LL>();
    omp_set_num_threads(n_threads);

    if(sequence_file == ""){
        cerr << "Error: Sequence file not given" << endl; exit(1);
    }

    if(out_prefix == ""){
        cerr << "Error: Prefix not given" << endl; exit(1);
    }

    if(baitfile == ""){
        cerr << "Error: Bait file not given" << endl; exit(1);
    }

    throwing_ofstream final_gaps_out(out_prefix + "-gaps.txt");
    throwing_ofstream final_crossings_out(out_prefix + "-crossings.txt");
    throwing_ofstream final_coverage_out(out_prefix + "-coverage.txt");
    throwing_ofstream cover_curve_out(out_prefix + "-cover-fractions.txt");
    throwing_ofstream match_positions(out_prefix + "-match-positions.txt");

    cerr << "Loading sequences" << endl;
    vector<string> sequences = read_sequences(sequence_file, true); // Reverse complements included to match the FM index. TODO: document.
    cerr << "Loading baits" << endl;
    vector<string> baits = read_sequences(baitfile, false); // Reverse complements not included

    assert(baits.size() > 0);
    for(LL i = 1; i < baits.size(); i++) assert(baits[i].size() == baits[i-1].size());
    
    LL bait_length = baits[0].size();

    FM_index fmi;
    if(fm_index_file == ""){
        cerr << "Constructing FM index" << endl;
        fmi.construct(sequences);
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

    callback_t write_coverage_callback = [&](const State& state){
        if(state.bait_id == baits.size()-1){
            cerr << "Writing coverage" << endl;
            for(auto& v : state.cover){
                for(LL b : v) final_coverage_out.stream << b << " ";
                final_coverage_out.stream << "\n";
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

    callback_t write_match_positions = [&](const State& state){
        // Line format: bait-id ref-id start-pos-in-ref mismatches is-rc-match
        for(pair<LL,LL> P : state.neighbors){
            LL doc_id, pos; std::tie(doc_id, pos) = P;
            LL dist_fw = hamming_distance_not_matching_N(state.bait, sequences[doc_id].substr(pos, bait_length));
            LL dist_rc = hamming_distance_not_matching_N(get_rc(state.bait), sequences[doc_id].substr(pos, bait_length));
            match_positions << state.bait_id << " " << doc_id << " " << pos << " " << min(dist_fw, dist_rc) << " " << (dist_rc < dist_fw) << "\n"; 
        }
    };

    vector<callback_t> callbacks = {gap_length_callback, write_coverage_callback, write_cover_curve_callback, count_branch_crossings_callback, write_match_positions};

    run(sequences, baits, fmi, d, g, bait_length, verbose, callbacks);


}