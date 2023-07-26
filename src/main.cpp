// Segment picture hybridization expansion
// Author: H. LaBollita (2023)

#include "configuration.hpp"
#include "segment.hpp"
#include "antisegment.hpp"
#include "hybridization.hpp"
#include "green.hpp"
#include "util.hpp"


int main(){

    double beta = 20;
    double h = 0;
    double t = 1.0;
    int    nt = 200;
    auto times = nda::zeros<double>(nt);
    for (int i=0; i < nt; i++){ times(i) = beta * (double)i/(nt-1); }

    Configuration c = Configuration(nda::vector<double>{1.0, 3.0},
                                    nda::vector<double>{2.0, 4.0}
            );

    GreensFunction g_ref = read_semi_circular_g_tau();
    Hybridization Delta = Hybridization(times, -0.25*t*t*g_ref.data, beta);

    auto A = nda::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0};
    std::cout << A << std::endl;
    auto B = roll(A, -2);
    std::cout << B << std::endl;
  return 0;
}
