#pragma once

#include "nda/nda.hpp"

template <typename T>
nda::vector<T> roll(nda::vector<T>& a, int shift) {
    auto N = a.size();
    auto A = nda::zeros<T>(N);
    for (int i=0; i < N; i++) {
        auto idx = i < N ? (N + i - shift) % N : i;
        A(idx) = a(i);
    }
    return A;
}

template <typename T>
nda::vector<T> hstack(nda::vector<T>& a, T& b) {
    auto N = a.size();
    auto A = nda::zeros<T>(N+1);
    for (auto i=0; i < N; i++) { A(i) = a(i); }
    A(N) = b;
    return A;
}

template <typename T>
nda::vector<T> stack(nda::vector<T>& a, nda::vector<T>& b) {
    auto N = a.size();
    auto M = b.size();
    auto A = nda::zeros<T>(N+M);
    for (auto i=0; i < N; i++) { A(i) = a(i); }
    for (auto j=N+1; j < M+N; j++) { A(j) = b(j); }
    return A;
}

int randomint(int min, int max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
}

