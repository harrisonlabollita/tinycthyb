#pragma once

#include "nda/nda.hpp"
#include <cmath>

template <typename T> nda::vector<T> deleteat(const nda::vector<T> old, int idx) {
  nda::vector<T> out = nda::zeros<T>(old.extent(0) - 1);
  int k = 0;
  for (int i = 0; i < old.size(); i++) {
    if (i != idx) {
      out(k) = old(i);
      k++;
    }
  }
  return out;
}

template <typename T> nda::vector<T> roll(const nda::vector<T> &a, int shift) {
  auto N = a.size();
  auto A = nda::zeros<T>(N);
  for (int i = 0; i < N; i++) {
    auto idx = i < N ? (N + i - shift) % N : i;
    A(idx) = a(i);
  }
  return A;
}

template <typename T> nda::vector<T> hstack(const nda::vector<T> &a, T &b) {
  auto N = a.size();
  auto A = nda::zeros<T>(N + 1);
  for (auto i = 0; i < N; i++) {
    A(i) = a(i);
  }
  A(N) = b;
  return A;
}

template <typename T>
nda::vector<T> stack(const nda::vector<T> &a, nda::vector<T> &b) {
  auto N = a.size();
  auto M = b.size();
  auto A = nda::zeros<T>(N + M);
  for (auto i = 0; i < N; i++) {
    A(i) = a(i);
  }
  for (auto j = N + 1; j < M + N; j++) {
    A(j) = b(j);
  }
  return A;
}

int randomint(int min, int max) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(min, max);
  return dis(gen);
}
template <typename T> int sign(T number) {
  return std::signbit(number) ? -1 : (number > 0 ? 1 : 0);
}
