#include "FM_index.hh"

bool operator==(const FM_index::Interval& A, const FM_index::Interval& B){
    return A.left == B.left && A.right == B.right;
}

std::ostream& operator<<(std::ostream& os, FM_index::Interval I){
    os << "[" << I.left << ", " << I.right << "]";
    return os;
}