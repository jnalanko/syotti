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
#include "FM_index.hh"
#include "cxxopts.hpp"

using namespace std;

int index_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Build an FM index.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("s,sequences", "Path to the fasta file of the input sequences.", cxxopts::value<string>()->default_value(""))
      ("o,output", "Path to the output FM-index file", cxxopts::value<string>()->default_value(""))
      ("h,help", "Print instructions.", cxxopts::value<bool>()->default_value("false"));

    cxxopts::ParseResult cli_params = options.parse(argc, argv);
    if (cli_params["help"].as<bool>() == true || original_argc == 1){
        cerr << options.help() << endl;
        return 1;
    }

    FM_index FMI;
    string infile = cli_params["sequences"].as<string>();
    string outfile = cli_params["output"].as<string>();
    throwing_ofstream out(outfile, ios_base::binary);
    vector<string> seqs = read_sequences(infile);
    FMI.construct(seqs);
    LL written = FMI.serialize(out.stream);
    cerr << "Wrote " << written << " bytes to " << outfile << endl;
}