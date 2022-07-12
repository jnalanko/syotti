#include <vector>
#include <string>
#include <algorithm>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include "input_reading.hh"
#include "cxxopts.hpp"

using namespace std;
typedef long long LL;

int subset_main(int argc, char** argv){
    cxxopts::Options options(argv[0], "Take a random subset of sequences in a fasta file.");
    int original_argc = argc; // It seems the CLI parsing library modifies argc, so store the original value

    options.add_options()
      ("i,input", "Input fasta file.", cxxopts::value<string>()->default_value(""))
      ("o,output", "Output fasta file.", cxxopts::value<string>()->default_value(""))
      ("n", "Size of the subset.", cxxopts::value<LL>()->default_value("-1"))
      ("h,help", "Print instructions.", cxxopts::value<bool>()->default_value("false"));

    cxxopts::ParseResult cli_params = options.parse(argc, argv);
    if (cli_params["help"].as<bool>() == true || original_argc == 1){
        cerr << options.help() << endl;
        return 1;
    }

    string in_filename = cli_params["input"].as<string>();
    string out_filename = cli_params["output"].as<string>();
    LL howmany = cli_params["n"].as<LL>();

    if(in_filename == ""){
        cerr << "Error: input filename not given" << endl;
        return 1;
    }

    if(out_filename == ""){
        cerr << "Error: output filename not given" << endl;
        return 1;
    }

    if(howmany == -1){
        cerr << "Error: subset size not given" << endl;
        return 1;
    }

    vector<pair<string,string> > seqs = parse_FASTA(in_filename);
    cerr << "Read " << seqs.size() << " sequences" << endl;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    cerr << "Randomizing" << endl;
    std::shuffle(seqs.begin(), seqs.end(), std::default_random_engine(seed));

    cerr << "Writing " << howmany << " sequences to " << out_filename << endl;
    throwing_ofstream out(out_filename);
    for(LL i = 0; i < howmany; i++){
        out << seqs[i].second << "\n" << seqs[i].first << "\n";
    }

    return 0;

}