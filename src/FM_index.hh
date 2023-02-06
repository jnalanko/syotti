#pragma once

#include <vector>
#include <cstdlib>
#include <string>
#include "sdsl/io.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "divsufsort.h"
#include "util.hh"

typedef int64_t LL;

using namespace std;

//extern "C" saint_t divsufsort(const sauchar_t *T, saidx_t *SA, saidx_t n);

// FM index for multiple input strings
// Is not aware for reverse compelements. If you want to index reverse compelements too,
// give that as input also
class FM_index{

public:

    struct Interval{
        LL left, right;
        LL size(){return right-left+1;}
        friend bool operator==(const Interval& A, const Interval& B);
        //friend std::ostream& operator<<(std::ostream& os, const Interval& I);
    };

private:

    sdsl::int_vector<40> SA; // Supports inputs up to 1TB
    sdsl::wt_hutu<sdsl::bit_vector> bwt;
    vector<int64_t> C_array;
    string concat;
    sdsl::bit_vector doc_starts;
    sdsl::rank_support_v<1> doc_starts_rs;
    sdsl::select_support_mcl<1> doc_starts_ss;
    Interval EMPTY_INTERVAL = {-1,0}; // Has length 0 by the formula end - start + 1
    Interval FULL_INTERVAL = {-1,0}; // Set after construction

    string concatenate(const vector<string>& inputs) const{
        string X;
        for(string S : inputs){
            for(char c : S) assert(c != '$' && c != '#');
            X += "$" + S;
        }
        X += "#";
        return X;
    }

    // No copying because of pointer members
    FM_index(const FM_index& other) = delete;
    FM_index& operator=(const FM_index& other) = delete;

public:

    Interval get_full_interval() const {return FULL_INTERVAL;};
    Interval get_empty_interval() const {return EMPTY_INTERVAL;};
    string get_concat() const {return concat;}
    LL SA_at(LL i) const {return SA[i];}
    LL C_array_at(char c) const {return C_array[c];}
    char bwt_at(LL i) const {return bwt[i];}
    LL bwt_rank(LL i, char c) const {return bwt.rank(i, c);}
    LL size() const {return concat.size();}

    FM_index(){}

    vector<int64_t> char_counts_to_C_array(const vector<int64_t>& counts) const{
        vector<int64_t> C(256); // Cumulative sum of counts

        // Compute cumulative sum of counts
        for(int64_t i = 0; i < (int64_t)C.size(); i++){
            C[i] = counts[i];
            if(i > 0) C[i] += C[i-1];
        }

        // Shift C to the right by one because that's how it's defined
        for(int64_t i = 256-1; i >= 0; i--){
            if(i == 0) C[i] = 0;
            else C[i] = C[i-1];
        }

        return C;

    }

    // Returns bytes written
    LL serialize(ostream& os){
        LL written = 0;

        LL C_array_size_bytes = C_array.size() * sizeof(int64_t);
        LL concat_size_bytes = concat.size() * sizeof(char);

        os.write((char*)&C_array_size_bytes, sizeof(C_array_size_bytes));
        os.write((char*)C_array.data(), C_array_size_bytes);

        os.write((char*)&concat_size_bytes, sizeof(concat_size_bytes));
        os.write((char*)concat.data(), concat_size_bytes);

        written += sizeof(C_array_size_bytes) + C_array_size_bytes;
        written += sizeof(concat_size_bytes) + concat_size_bytes;

        written += SA.serialize(os);
        written += bwt.serialize(os);
        written += doc_starts.serialize(os);
        written += doc_starts_rs.serialize(os);
        written += doc_starts_ss.serialize(os);
        return written;
    }

    void load(istream& is){
        LL C_array_size_bytes, concat_size_bytes;
        
        // Read C array
        is.read((char*)&C_array_size_bytes, sizeof(C_array_size_bytes));
        C_array.resize(C_array_size_bytes / sizeof(int64_t));
        is.read((char*)C_array.data(), C_array_size_bytes);

        // Read concat
        is.read((char*)&concat_size_bytes, sizeof(concat_size_bytes));
        concat.resize(concat_size_bytes / sizeof(char), '\0');
        is.read((char*)concat.data(), concat_size_bytes);

        // Read sdsl structures
        SA.load(is);
        bwt.load(is);
        doc_starts.load(is);
        doc_starts_rs.load(is);
        doc_starts_ss.load(is);

        // Set supports
        doc_starts_rs.set_vector(&doc_starts);
        doc_starts_ss.set_vector(&doc_starts);

        // Set full interval
        FULL_INTERVAL = {(LL)0, (LL)bwt.size()-1};
    }

    void construct(const vector<string>& inputs){
        concat = concatenate(inputs);

        this->doc_starts.resize(concat.size());
        for(int64_t i = 0; i < concat.size(); i++){
            doc_starts[i] = (concat[i] == '$' ? 1 : 0);
        }
        sdsl::util::init_support(doc_starts_rs, &doc_starts);
        sdsl::util::init_support(doc_starts_ss, &doc_starts);

        // Compute the suffix array
        int64_t* temp_SA = (int64_t*)malloc(sizeof(int64_t) * concat.size());
        divsufsort64((sauchar_t*)(&concat[0]), temp_SA, concat.size());
        SA.resize(concat.size());
        for(LL i = 0; i < concat.size(); i++) SA[i] = temp_SA[i];
        free(temp_SA);

        string BWT_string(concat.size(), '.');
        for(int64_t i = 0; i < concat.size(); i++){
            if(SA[i] == 0) BWT_string[i] = concat.back();
            else BWT_string[i] = concat[SA[i]-1];
        }
        construct_im(this->bwt, (const char*)(&BWT_string[0]), 1);

        vector<int64_t> char_counts(256);
        for(char c : BWT_string) char_counts[c]++;
        C_array = char_counts_to_C_array(char_counts);
        
        FULL_INTERVAL = {(LL)0, (LL)bwt.size()-1};
    }

    // Left extend the input interval with c
    Interval left_extend(Interval I, char c) const{
        if(I == EMPTY_INTERVAL) return I;
        int64_t r1 = bwt.rank(I.left, c);
        int64_t r2 = bwt.rank(I.right + 1,c);
        Interval I_new = {C_array[c] + r1, C_array[c] + r2 - 1};
        if(I_new.left > I_new.right) return EMPTY_INTERVAL;
        else return I_new;
    }

    Interval search(const string& x) const{
        Interval I = FULL_INTERVAL;
        for(LL i = (LL)x.size()-1; i >= 0; i--){
            I = left_extend(I, x[i]);
            if(I == EMPTY_INTERVAL) return EMPTY_INTERVAL;
        }
        return I;
    }

    // Returns (doc id, position in doc) of the suffix starting at lex_rank
    // If the suffix starts with a dollar, returns (-1,-1)
    pair<int64_t, int64_t> locate(LL lex_rank) const{
        LL p = SA[lex_rank];
        if(concat[p] == '$') return {-1,-1};
        int64_t doc_id = doc_starts_rs(p) - 1;
        int64_t offset = p - (doc_starts_ss(doc_id+1) + 1);
        return {doc_id, offset};
    }

    // Returns pairs (doc_id, pos) such that if X is the document with id doc_id then X[pos..pos+L-1] = S.
    vector<pair<int64_t, int64_t> > subpattern_search(const string& S, LL L, LL pos) const{
        assert(L > 0 && S.size() <= L);
        vector<pair<int64_t, int64_t> > ans;
        Interval I = FULL_INTERVAL;
        for(LL i = (LL)S.size()-1; i >= 0; i--){
            I = left_extend(I, S[i]);
            if(I == EMPTY_INTERVAL) return ans;
        }

        // I is now the lexicographic range of S
        for(LL i = I.left; i <= I.right; i++){
            LL p = (LL)SA[i] - pos;
            if(p < 0) continue;
            string candidate = concat.substr(p,L);
            if(candidate.size() < L) continue;
            bool good = true;
            for(char c : candidate) if(c == '$' || c == '#') good = false;
            if(good){
                int64_t doc_id = doc_starts_rs(p) - 1;
                int64_t offset = p - (doc_starts_ss(doc_id+1) + 1);
                ans.push_back({doc_id, offset});
            }
        }
        return ans;
    }

    bool operator==(const FM_index& other) const{
        const FM_index& fmi1 = *this;
        const FM_index& fmi2 = other;

        bool same = true;

        same &= fmi1.SA == fmi2.SA;
        same &= fmi1.C_array == fmi2.C_array;
        same &= fmi1.concat == fmi2.concat;
        same &= fmi1.doc_starts == fmi2.doc_starts;

        same &= fmi1.bwt.size() == fmi2.bwt.size();
        for(LL i = 0; i < fmi1.bwt.size(); i++)
            same &= fmi1.bwt[i] == fmi2.bwt[i];

        for(LL i = 0; i <= fmi1.doc_starts.size(); i++)
            same &= fmi1.doc_starts_rs(i) == fmi2.doc_starts_rs(i);

        for(LL i = 1; i <= fmi1.doc_starts_rs.rank(fmi1.doc_starts.size()); i++)
            same &= fmi1.doc_starts_ss(i) == fmi2.doc_starts_ss(i);

        same &= (fmi1.FULL_INTERVAL == fmi2.FULL_INTERVAL);
        return same;

    }

};

bool operator==(const FM_index::Interval& A, const FM_index::Interval& B);
std::ostream& operator<<(std::ostream& os, FM_index::Interval I);