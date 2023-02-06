
#include "util.hh"

char get_rc(char c){
    if(c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N'){
        cerr << "Warning: non-ACGTN character found: " << c << endl;
    }
    if(c == 'A') return 'T';
    else if(c == 'C') return 'G';
    else if(c == 'G') return 'C';
    else if(c == 'T') return 'A';
    else return 'N'; // The reverse complement of N is N
}

string get_rc(string S){
    std::reverse(S.begin(), S.end());
    for(char& c : S){
        c = get_rc(c);
    }
    return S;
}

string get_canonical(const string& S){
    string T = get_rc(S);
    return S < T ? S : T;
}


vector<string> split(string s, char delimiter){
    stringstream test(s);
    string segment;
    vector<string> seglist;

    while(getline(test, segment, delimiter)){
        seglist.push_back(segment);
    }
    return seglist;
}

// Get smallest g-mer of S
string get_minimizer(const string& S, LL g){
    assert(S.size() >= g);
    string M = get_canonical(S.substr(0,g));
    for(LL i = 1; i < S.size()-g+1; i++){
        string x = get_canonical(S.substr(i,g));
        if(x < M) M = x;
        
    }
    return M;
}

// 'N' is defined to not match to anything
LL hamming_distance_not_matching_N(const string& A, const string& B){
    assert(A.size() == B.size());
    LL ans = 0;
    for(LL i = 0; i < A.size(); i++){
        if(A[i] == 'N' || B[i] == 'N' || A[i] != B[i]) ans++;
    }
    return ans;
}

vector<string> read_sequences(string filename, bool append_reverse_complements){
    vector<string> seqs;
    SeqIO::Reader<> reader(filename);
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        seqs.push_back(string(reader.read_buf));
    }

    if(append_reverse_complements){
        LL n = seqs.size();
        for(LL i = 0; i < n; i++)
            seqs.push_back(get_rc(seqs[i]));
    }
    return seqs;
}

vector<string> get_all_distinct_kmers_including_rc(const vector<string>& v, LL k){
    cerr << "Listing distinct k-mers (rc included)" << endl;
    Progress_printer pp(v.size(), 100);
    unordered_set<string> hashset;
    for(string S : v){
        for(LL i = 0; i < S.size()-k+1; i++){
            hashset.insert(S.substr(i,k));
            hashset.insert(get_rc(S.substr(i,k)));
        }
        pp.job_done(to_string(k) + "-mers");
    }
    vector<string> ans;
    for(string S : hashset) ans.push_back(S);
    return ans;
}

void check_writable(string filename){
    ofstream F(filename);
    if(!F.good()){
        cerr << "Error writing to file: " << filename << endl;
        exit(1);
    }
}
