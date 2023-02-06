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

TEST(greedy_test, small_hand_crafted){
    vector<string> seqs = {"AAAAAACCCCCCATATATAGTTTTTTTT",
                           "NNNNNNNNNNNNN",
                           "AAAAAAAACTATATATGGGGGGTTTTTT", // First again
                           "AAAAAAAACTATATATGGGGGGTTTTTT", // RC of the first
                           "AANNAANNCCNNATNNATNNTTNNTTNN", // The first but N's such that there is no common 5-mer with 1 mismatch
                           "AANNAANNCTNNATNNGGNNGGNNTTNN", // The RC of the first but N's such that there is no common 5-mer with 1 mismatch
                           "NNNNNATATATANNNNNNNNN", // Island in the middle should be covered by bait TATAT
                           "TACGT", // Unique
                           "ACGTA", // RC of above
                           "TGTXT", // Non-ACGTN character
                           "AYACA", // RC of above with X changed to Y
                           };
    LL d = 1;
    LL g = 2;
    LL bait_length = 5;

    FM_index fmi;
    fmi.construct(seqs);
    FM_NeighborCandidateFunction FM_NCF;
    FM_NCF.init(&fmi, g);
    NeighborFunction FM_NF;
    FM_NF.init(&FM_NCF, &seqs, d, bait_length);

    Greedy G_FM;
    G_FM.init(&FM_NF, &seqs, bait_length, d, g, false, 1);
    Greedy::Result result = G_FM.run();

    // Covering the first sequence should happen like this:
    // AAAAA covers the prefix AAAAAAC and by reverse complement the suffix GTTTTTTTT.
    // CCCCC covers CCCCCCA.
    // TATAT covers TATATA (also the last A because of reverse complement ATATA).
    // The Ns don't match to each other so they are all covered separately.
    // The third input sequence is just the reverse complement of the first, so it is automatically covered.
    vector<string> expected_baits = {"AAAAA","CCCCC","TATAT", // First sequence
                                     "NNNNN","NNNNN","NNNNN", // Second sequence
                                                              // 3. sequence: already covered
                                                              // 4. sequence: already covered
                                     "AANNA","ANNCC","NNATN","NATNN","TTNNT","NTTNN",  // 5. sequence
                                     "AANNA","ANNCT","NNATN","NGGNN","GGNNT","NTTNN", // 6. sequence
                                     "NNNNN","NNNNN","NNNNN", // 7. sequence
                                     "TACGT", // 8. sequence
                                              // 9. sequence: already covered as RC of 8.
                                     "TGTNT", // 10. sequence. The x should be replaced with an N
                                     "", // 10. sequence: already covered as RC of 9
                                      };

    cout << result.baits << endl;
    cout << expected_baits << endl;
    ASSERT_TRUE(result.baits == expected_baits);
}

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