// Segment picture hybridization expansion
// Author: H. LaBollita (2023)

#include "configuration.hpp"
#include "segment.hpp"
#include "antisegment.hpp"
#include "hybridization.hpp"

int main(){

    double beta = 20;
    double h = 0;
    double t = 1.0;
    int    nt = 200;
    auto times = nda::zeros<double>(nt);
    for (int i=0; i < nt; i++){
        times(i) = beta * (double)i/(nt-1) ;
    }

    auto g_ref = nda::zeros<double>(nt); 

    Configuration c = Configuration(nda::vector<double>{1.0, 3.0},
                                    nda::vector<double>{2.0, 4.0}
            );

    Hybridization Delta = Hybridization(times, g_ref, beta);
    auto d = Delta(times);

  return 0;
}
