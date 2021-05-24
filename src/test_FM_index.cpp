#include "FM_index.hh"
#include "stdlib_printing.hh"
#include "util.hh"
#include "throwing_streams.hh"
#include <gtest/gtest.h>
#include <cassert>

TEST(FM_index_test, serialization){
    vector<string> seqs;
    seqs.push_back("ABBABCABACBACB");
    seqs.push_back("AABAB");
    seqs.push_back("BACBACBA");

    FM_index fmi1;
    fmi1.construct(seqs);
    throwing_ofstream out("temp/test.fmi", ios_base::binary);
    fmi1.serialize(out.stream);
    out.close();

    FM_index fmi2;
    throwing_ifstream in("temp/test.fmi", ios_base::binary);
    fmi2.load(in.stream);

    for(LL i = 0; i < min((LL)10, (LL)(fmi1.size())); i++){
        cout << fmi1.SA_at(i) << " " << fmi2.SA_at(i) <<  endl;
    }

    ASSERT_TRUE(fmi1 == fmi2);

    cerr << "Serialization ok" << endl;
}

TEST(FM_index_test, search){
    FM_index FMI;
    vector<string> test;
    test.push_back("ABBABCABACBACB");
    test.push_back("AABAB");
    test.push_back("BACBACBA");
    FMI.construct(test);

    string concat = FMI.get_concat();
    cout << concat << endl;

    // Check that the suffix array is correct
    for(LL i = 0; i < concat.size(); i++){
        cout << i << " " << concat.substr(FMI.SA_at(i)) << endl;
        if(i > 0) ASSERT_TRUE(concat.substr(FMI.SA_at(i)) > concat.substr(FMI.SA_at(i-1)));
    }

    // Test searching for XAACB
    FM_index::Interval I = FMI.get_full_interval();
    FM_index::Interval I_B = FMI.left_extend(I,'B');
    FM_index::Interval I_CB = FMI.left_extend(I_B,'C');
    FM_index::Interval I_ACB = FMI.left_extend(I_CB,'A');
    FM_index::Interval I_AACB = FMI.left_extend(I_ACB, 'A');
    FM_index::Interval I_XAACB = FMI.left_extend(I_AACB, 'X');
    cout << I_B << " " << I_CB << " " << I_ACB << " " << I_AACB << " " << I_XAACB << endl;
    ASSERT_TRUE((I_B == (FM_index::Interval){15,25}));
    ASSERT_TRUE((I_CB == (FM_index::Interval){27,30}));
    ASSERT_TRUE((I_ACB == (FM_index::Interval){11,14}));
    ASSERT_TRUE((I_AACB == FMI.get_empty_interval()));
    ASSERT_TRUE((I_XAACB == FMI.get_empty_interval()));

    // Test subpattern search
    vector<pair<LL,LL> > hits = FMI.subpattern_search("BAC",6,2);
    for(pair<LL,LL> H : hits) {
        string S = test[H.first].substr(H.second, 6);
        cout << S << endl;
        ASSERT_TRUE(S.substr(2,3) == "BAC");
    }
    cout << FMI.subpattern_search("ACBA", 10, 2) << endl;
    ASSERT_TRUE(FMI.subpattern_search("ACBA", 10, 2).size() == 0); // ACBA exists but there will be dollars

    // Test locating
    FM_index::Interval I_CABACBAC = FMI.search("CABACBAC");
    cout << I_CABACBAC << endl;
    ASSERT_TRUE(I_CABACBAC.left == I_CABACBAC.right);
    ASSERT_TRUE(I_CABACBAC.left == 26);
    LL doc_id, pos;
    std::tie(doc_id, pos) = FMI.locate(I_CABACBAC.left);
    cout << doc_id << " " << pos << endl;
    ASSERT_TRUE(doc_id == 0);
    ASSERT_TRUE(pos == 5);
}