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
#include <omp.h>
#include "progress.hh"
#include "util.hh"
#include "greedy.hh"
#include "cxxopts.hpp"

using namespace std;

int main(int argc, char** argv){

    cxxopts::Options options(argv[0],"Computes a greedy bait cover.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("L,bait-len", "Length of the baits.", cxxopts::value<LL>()->default_value("120"))
      ("d,hamming-distance", "Number of allowed mismatches in the baits.", cxxopts::value<LL>()->default_value("40"))
      ("s,sequences", "Path to a fasta file of the input sequences.", cxxopts::value<string>()->default_value(""))
      ("fm-index-out", "The algorithm is based on FM-index, which we build at the start. Building the index can take a lot of time and memory. Use this option to save the FM-index to disk so that you can later run the algorithm with different parameters re-using the same FM-index. (optional).", cxxopts::value<string>()->default_value(""))
      ("f,fm-index", "Path to a previously saved FM-index on disk (--fm-index-out). This option loads the FM index from disk instead of building it again.", cxxopts::value<string>()->default_value(""))
      ("o,out", "Filename prefix for the output files.", cxxopts::value<string>()->default_value(""))
      ("r,randomize", "Randomize the processing order of the sequences in the greedy algorithm.", cxxopts::value<bool>()->default_value("false"))
      ("t,n-threads", "Maximum number of parallel threads. The program is not very well optimized for parallel processing, so don't expect much of a speedup here.", cxxopts::value<LL>()->default_value("1"))
      ("c,cutoff", "Stop the greedy algorithm after this fraction of positions is covered. For example: 0.99.", cxxopts::value<double>()->default_value("1"))
      ("g,seed-len", "The length of the seeds in the FM-index seed-and-extend approximate string search subroutine. A lower value will find more matches, but will be slower.", cxxopts::value<LL>()->default_value("20"))
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
    double cutoff = cli_params["cutoff"].as<double>();
    bool randomize = cli_params["randomize"].as<bool>();

    omp_set_num_threads(n_threads);

    string out_prefix = cli_params["out"].as<string>();
    string baits_outfile = out_prefix + "-baits.fna";
    string cover_fractions_outfile = out_prefix + "-cover-fractions.txt";
    string cover_marks_outfile = out_prefix + "-cover-marks.txt";
    check_writable(baits_outfile);
    check_writable(cover_fractions_outfile);
    check_writable(cover_marks_outfile);

    vector<string> seqs = read_sequences(cli_params["sequences"].as<string>(), true); // Appends reverse complement
    for(string& seq : seqs){
        if(seq.size() < bait_length){
            cerr << "Error: sequence is shorter than bait length " << bait_length << ": " << seq << endl;
            return 1;
        }
    }

    string fmi_path = cli_params["fm-index"].as<string>();
    string fmi_out = cli_params["fm-index-out"].as<string>();
    FM_index fmi;
    if(fmi_path == ""){
        cerr << "Constructing FM index" << endl;
        fmi.construct(seqs);
        cerr << "FM index construction done" << endl;
        if(fmi_out != ""){
            cerr << "Saving the FM-index to " << fmi_out << endl;
            throwing_ofstream fmi_out_stream(fmi_out, ios_base::binary);
            fmi.serialize(fmi_out_stream.stream);
        }
    } else{
        throwing_ifstream in(fmi_path, ios_base::binary); // Throws on error
        cerr << "Loading FM-index" << endl;
        fmi.load(in.stream);
        cerr << "FM index loaded" << endl;
    }

    FM_NeighborCandidateFunction NCF;
    NCF.init(&fmi, g);
    NeighborFunction NF;
    NF.init(&NCF, &seqs, d, bait_length);
    Greedy alg;
    alg.init(&NF, &seqs, bait_length, d, g, randomize, cutoff);
    Greedy::Result res = alg.run();

    cerr << res.baits.size() << " baits" << endl;
    
    throwing_ofstream baits_out(baits_outfile);
    throwing_ofstream cover_marks_out(cover_marks_outfile);
    throwing_ofstream cover_fractions_out(cover_fractions_outfile);

    cerr << "Bait output file: " << baits_outfile << endl;
    cerr << "Cover fractions output file: " << cover_fractions_outfile << endl;
    cerr << "Cover marks output file: " << cover_marks_outfile << endl;

    for(string& S : res.baits) baits_out << ">\n" << S << "\n";
    for(double x : res.cover_fractions) cover_fractions_out << x << "\n";
    for(vector<bool>& v : res.covered){
        for(bool b : v) cover_marks_out << b;
        cover_marks_out << "\n";
    }
    
}
