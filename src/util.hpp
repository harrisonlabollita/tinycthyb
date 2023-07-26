#pragma once

#include "nda/nda.hpp"

template <typename T>
nda::vector<T> roll(nda::vector<T>& a, int shift) {
    auto N = a.size();
    auto A = nda::zeros<double>(N);
    for (int i=0; i < N; i++) {
        auto idx = i < N ? (N + i - shift) % N : i;
        A(idx) = a(i);
    }
    return A;
}
