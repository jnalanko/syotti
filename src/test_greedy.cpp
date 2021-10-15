#include "greedy.hh"
#include "stdlib_printing.hh"
#include <cassert>
#include <gtest/gtest.h>

TEST(greedy_test, reverse_complement){
    string S = "ATGNAC";
    string rev_S = "GTNCAT"; // Reverse complement of N is N
    ASSERT_TRUE(get_rc(S) == rev_S);
}


TEST(greedy_test, hamming_distance_not_matching_N){
    string S = "AACCGGTTNN";
    string T = "ATCTGTTANA";
    ASSERT_TRUE(hamming_distance_not_matching_N(S,T) == 6);
}

/*
TEST(greedy_test, small_hand_crafted){
    vector<string> seqs = {"AAAAAACCCCCCAAAA","NNNNNNNNNNNNN"};
    LL d = 1;
    LL g = 3;
    LL bait_length = 5;

    FM_index fmi;
    fmi.construct(seqs);
    FM_NeighborCandidateFunction FM_NCF;
    FM_NCF.init(&fmi, g);
    NeighborFunction FM_NF;
    FM_NF.init(&FM_NCF, &seqs, d, bait_length);

    Greedy G_FM;
    G_FM.init(&FM_NF, &seqs, bait_length, d, g, false, 1);

    // TODO

}*/

TEST(greedy_test, FM_against_hash){
    LL d = 4;
    LL g = 12;
    LL bait_length = 20;
    vector<string> seqs = read_sequences("testcases/coli3.fna");
    FM_index fmi;
    fmi.construct(seqs);
    FM_NeighborCandidateFunction FM_NCF;
    FM_NCF.init(&fmi, g);
    NeighborFunction FM_NF;
    FM_NF.init(&FM_NCF, &seqs, d, bait_length);

    HashNeighborCandidateFunction Hash_NCF;
    Hash_NCF.init(seqs, bait_length, g);
    NeighborFunction Hash_NF;
    Hash_NF.init(&Hash_NCF, &seqs, d, bait_length);

    Greedy G_FM;
    G_FM.init(&FM_NF, &seqs, bait_length, d, g, false, 1);

    Greedy G_hash;
    G_hash.init(&Hash_NF, &seqs, bait_length, d, g, false, 1);

    Greedy::Result res_hash = G_hash.run();
    Greedy::Result res_FM = G_FM.run();
    cout << res_hash.baits.size() << " " << res_FM.baits.size() << endl;
    ASSERT_TRUE(res_hash.baits.size() == res_FM.baits.size());
}