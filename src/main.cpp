// Segment picture hybridization expansion
// Author: H. LaBollita (2023)

#include "configuration.hpp"
#include "segment.hpp"
#include "antisegment.hpp"
#include "hybridization.hpp"
#include "green.hpp"
#include "solver.hpp"


int main(){

    double beta = 20;
    double h = 0;
    int    nt = 200;
    auto times = nda::zeros<double>(nt);
    for (int i=0; i < nt; i++){ times(i) = beta * (double)i/(nt-1); }

    GreensFunction g_ref = read_semi_circular_g_tau();
    Hybridization Delta = Hybridization(times, -0.25*g_ref.data, beta);
    auto e = Expansion(beta, h, Delta);
    auto c = Configuration(nda::vector<double>{1.0}, nda::vector<double>{3.0} );
    auto d = Determinant(c, e);
    auto t = trace(c, e);

    std::cout << "d.mat = " << d.mat << std::endl;
    std::cout << "d.val = " << d.value << std::endl;
    std::cout << "trace = " << t << std::endl;

  return 0;
}
