#include <vector>
#include <string>
#include <algorithm>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include "input_reading.hh"

using namespace std;
typedef long long LL;

int main(int argc, char** argv){
    if(argc == 1){
        cerr << "Usage in.fna out.fna howmany" << endl;
        return 1;
    }
    string in_filename = argv[1];
    string out_filename = argv[2];
    LL howmany = stoll(argv[3]);
    vector<pair<string,string> > seqs = parse_FASTA(in_filename);
    cerr << "Shuffling " << seqs.size() << " sequences" << endl;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(seqs.begin(), seqs.end(), std::default_random_engine(seed));

    cerr << "Writing the first " << howmany << " to " << out_filename << endl;
    throwing_ofstream out(out_filename);
    for(LL i = 0; i < howmany; i++){
        out << seqs[i].second << "\n" << seqs[i].first << "\n";
    }

}