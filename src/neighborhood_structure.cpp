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
#include "other_greedy.hh"
#include "cxxopts.hpp"
#include "neighborhood_structure.hh"

using namespace std;

int main(int argc, char** argv){

    cxxopts::Options options("Compute neighborhood structure","A precomputed FM-index can be given.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("f,fm-index", "Path to the FM index file of the sequences (optional)", cxxopts::value<string>()->default_value(""))
      ("fm-index-out", "Path where to save the FM-index (optional)", cxxopts::value<string>()->default_value(""))
      ("s,sequences", "Path to a fasta file of the input sequences", cxxopts::value<string>()->default_value(""))
      ("o,out", "Filename of the outputfile)", cxxopts::value<string>()->default_value(""))
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
    string outfile = cli_params["out"].as<string>();
    check_writable(outfile);

    vector<string> seqs = read_sequences(cli_params["sequences"].as<string>(), true); // Also appends reverse complements
    cerr << "Read " << seqs.size() << " sequences (reverse complements included)" << endl;

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

     Neighborhood_Structure NS;
     Neighborhood_Structure_Builder NSB;
     throwing_ofstream os(outfile, ios_base::binary);

     NSB.build(NS, &fmi, &seqs, bait_length, d, g);
     NS.serialize(os.stream);


}