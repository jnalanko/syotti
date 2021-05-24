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

using namespace std;

int main(int argc, char** argv){
    if(argc != 3){
        cerr << "Constructs the FM index of the input sequences and their reverse complements" << endl;
        cerr << "See FM_index.hh for details" << endl;
        cerr << "Appends reverse complement" << endl; // Todo: update this is changes
        cerr << "Usage: " << argv[0] << " data.fna out.fmi" << endl;
        return 1;
    }

    FM_index FMI;
    string outfile = argv[2];
    throwing_ofstream out(outfile, ios_base::binary);
    vector<string> seqs = read_sequences(argv[1], true); // Appends reverse complements also
    FMI.construct(seqs);
    LL written = FMI.serialize(out.stream);
    cerr << "Wrote " << written << " bytes to " << outfile << endl;
}