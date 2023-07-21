#pragma once

#include <algorithm>
#include <nda/nda.hpp>

class Hybridization {
    private:
        nda::vector<double> times;
        nda::vector<double> values;
        double beta;
    public:
        Hybridization(const nda::vector<double>& times,
                      const nda::vector<double>& values,
                      double beta) : times(times), values(values), beta(beta) {}

        double operator()(double time) const {

            auto call = [this](double t) 


            double s = 1.0;
            if (time < 0.0)  {
                s = -1.0;
                time += beta;
            }

        }

}
