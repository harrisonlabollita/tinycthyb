#pragma once

#include <cmath>
#include <fstream>
#include <nda/nda.hpp>
#include <string>

class GreensFunction {

public:
  double beta;
  const nda::vector<double> data;
  double sign;
  int N;

  GreensFunction(double beta, int N)
      : beta(beta), data(nda::zeros<double>(N)), sign(0.0), N(N) {}

  GreensFunction(double beta, const nda::vector<double> &data)
      : beta(beta), data(data), N(data.shape()[0]), sign(0.0) {}

  GreensFunction(double beta, const nda::vector<double> &data, double sign)
      : beta(beta), data(data), sign(sign), N(data.shape()[0]) {}

  int length() { return N; }

  GreensFunction operator*(double scalar) {
    auto new_data = data * scalar;
    return GreensFunction(beta, new_data, sign);
  }

  GreensFunction operator/(double scalar) {
    auto new_data = data / scalar;
    return GreensFunction(beta, new_data, sign);
  }

  void accumulate(double time, double value) {
    if (time < 0.0) {
      value *= -1;
      time += beta;
    }
    auto idx = static_cast<int>(std::floor(length() * time / beta));
    data(idx) += value;
  }

  int write_data(std::string filename) const {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
      std::cerr << "Failed to open file!" << std::endl;
      return 1;
    }

    double dt = beta / N;

    for (auto val : data) {
      outputFile << val / (-sign * beta * dt) << " ";
    }
    outputFile.close();
    return 0;
  }
};

GreensFunction read_semi_circular_g_tau(void) {

  double beta = 20;
  int nt = 200;
  auto data = nda::zeros<double>(nt);
  double sign = 0.0;

  std::fstream inputFile("gref.txt");

  std::string line;
  int idx = 0;
  while (std::getline(inputFile, line)) {
    double value = std::stod(line);
    data(idx) = value;
    idx += 1;
  }
  inputFile.close();

  return GreensFunction(beta, data, sign);
}
