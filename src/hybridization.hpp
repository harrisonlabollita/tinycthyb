#pragma once

#include <algorithm>
#include <nda/nda.hpp>

namespace tinycthyb {

class Hybridization {
public:
  nda::vector<double> times;
  nda::vector<double> values;
  double beta;

  Hybridization(nda::vector<double> times, nda::vector<double> values,
                double beta)
      : times(times), values(values), beta(beta) {}

  double operator()(double t) {
    double s = 1.0;
    if (t < 0.0) {
      s = -1.0;
      t += beta;
    }

    auto it = std::lower_bound(times.begin(), times.end(), t); // iterator to element
    int idx = std::distance(times.begin(), it);
    idx = idx == 0 ? 1 : idx;
    double ti = times(idx - 1);
    double tf = times(idx);
    double vi = values(idx - 1);
    double vf = values(idx);
    return s * (vi + (t - ti) * (vf - vi) / (tf - ti));
  }

  nda::vector<double> operator()(nda::vector<double> time) {
    auto out = nda::zeros<double>(time.size());
    for (int i = 0; i < time.size(); i++) {
      out(i) = (*this)(time(i));
    }
    return out;
  }
};

} // namespace tinycthyb
