#pragma once

#include "iostream"
#include <nda/nda.hpp>


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
    int length() const { return t_i.shape()[0]; }
    void print() const { std::cout << "Configuration(" << t_i << ", " << t_f << ")" << std::endl; }
};

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


