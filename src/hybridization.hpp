#pragma once

#include <algorithm>
#include <nda/nda.hpp>

class Hybridization {
    public:
        nda::vector<double> times;
        nda::vector<double> values;
        double beta;

        Hybridization(nda::vector<double> times,
                      nda::vector<double> values,
                      double beta) : 
            times(times), values(values), beta(beta) {}

        double operator()(double time) const {

            double s = 1.0;
            if (time < 0.0) {
                s = -1.0;
                time += beta;
            }

            auto it = std::lower_bound(times.begin(), times.end(), time); // iterator to element
            int idx = std::distance(times.begin(), it);
            idx = idx == 0 ? 1 : idx;
            std::cout << "idx = " << idx << std::endl;
            double ti = times(idx-1);
            double tf = times(idx);
            double vi = values(idx-1);
            double vf = values(idx);
            return s * (vi + (time-ti)*(vf-vi) /(tf-ti));
            
        }

        nda::vector<double> operator()(const nda::vector<double> times){
            auto out = nda::zeros<double>(times.shape()[0]);
            for(int i=0; i < times.extent(0); i++) {out(i) = (*this)(times[i]); }
            return out;
       }
};
