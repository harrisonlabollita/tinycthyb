#pragma once

#include "iostream"
#include <nda/nda.hpp>

#include "util.hpp"

template <typename T>
nda::vector<T> deleteat(nda::vector<T> old, int idx) {
    nda::vector<T> out = nda::zeros<T>(old.extent(0)-1);
    int k=0;
    for (int i =0; i < old.extent(0); i++) {
        if (i != idx) { 
            out[k] = old[i]; 
            k++;
        }

    }
    return out;
}


struct InsertMove {
    double t_i;
    double t_f;
    double l;
    InsertMove(double t_i, double t_f, double l) 
        : t_i(t_i), t_f(t_f), l(l) {}
};

struct RemovalMove {
    int i_idx;
    int f_idx;
    double l;
    RemovalMove(int i_idx, int f_idx, double l) 
        : i_idx(i_idx), f_idx(f_idx), l(l) {}
};

struct Configuration {
    nda::vector<double> t_i; 
    nda::vector<double> t_f;

    Configuration(nda::vector<double> i, nda::vector<double> f){
        if (i.shape()[0] == 0) {
            t_i = nda::vector<double>{};
            t_f = nda::vector<double>{};
        }
        std::sort(i.begin(), i.end());
        std::sort(f.begin(), f.end());
        t_i = i;
        t_f = f;
    }

    // Assignment operator for copying
    Configuration& operator=(const Configuration& other) {
        if (this != &other) {
            t_i = other.t_i;
            t_f = other.t_f;
        }
        return *this;
    }

    int length() const { return t_i.size(); }
    void print() const { std::cout << "Configuration(" << t_i << ", " << t_f << ")" << std::endl; }


    Configuration operator+(InsertMove& move) {
        nda::vector<double> new_t_i = t_i;
        nda::vector<double> new_t_f = t_f;
        new_t_i = hstack(new_t_i, move.t_i);
        new_t_f = hstack(new_t_f, move.t_f);
        return Configuration(new_t_i, new_t_f);
    }

    Configuration operator+(RemovalMove& move) {
        if (move.i_idx >= 0 && length() > 0) {
            nda::vector<double> new_t_i = t_i;
            nda::vector<double> new_t_f = t_f;
            new_t_i = deleteat(new_t_i, move.i_idx);
            new_t_f = deleteat(new_t_f, move.f_idx);
            return Configuration(new_t_i, new_t_f);
        } else {
            return *this;
        }
    }
};

