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
