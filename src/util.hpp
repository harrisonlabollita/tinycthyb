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

int randomint(int min, int max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
}

