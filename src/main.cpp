// Segment picture hybridization expansion
// Author: H. LaBollita (2023)

#include "configuration.hpp"
#include "segment.hpp"
#include "antisegment.hpp"
#include "hybridization.hpp"
#include "green.hpp"
#include "solver.hpp"
#include "util.hpp"


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
    std::cout << "is_segment_proper = " << is_segment_proper(c) << std::endl;

    c = Configuration(nda::vector<double>{19.0}, nda::vector<double>{1.0} );
    d = Determinant(c, e);
    t = trace(c, e);

    std::cout << "d.mat = " << d.mat << std::endl;
    std::cout << "d.val = " << d.value << std::endl;
    std::cout << "trace = " << t << std::endl;
    std::cout << "is_segment_proper = " << is_segment_proper(c) << std::endl;

    c = Configuration(nda::vector<double>{1.0, 5.0}, nda::vector<double>{3.0, 7.0} );
    d = Determinant(c, e);
    t = trace(c, e);

    std::cout << "d.mat = " << d.mat << std::endl;
    std::cout << "d.val = " << d.value << std::endl;
    std::cout << "trace = " << t << std::endl;
    std::cout << "is_segment_proper = " << is_segment_proper(c) << std::endl;

    c = Configuration(nda::vector<double>{2.0, 19.0}, nda::vector<double>{5.0, 1.0} );
    d = Determinant(c, e);
    t = trace(c, e);

    std::cout << "d.mat = " << d.mat << std::endl;
    std::cout << "d.val = " << d.value << std::endl;
    std::cout << "trace = " << t << std::endl;
    std::cout << "is_segment_proper = " << is_segment_proper(c) << std::endl;

    c = Configuration(nda::vector<double>{1.0, 3.0, 5.0}, nda::vector<double>{0.5, 2.0, 4.0} );
    std::cout << "is_segment_proper = " << is_segment_proper(c) << std::endl;
    
    std::cout << "Segmnets" << std::endl;
    for (auto s : segments(c) ) {
        s.print();
    }

    std::cout << "Anti-segmnets" << std::endl;
    for (auto s : antisegments(c) ) {
        s.print();
    }

    auto s = onsegment(1.5, c);
    if (s) { (*s).print(); 
    } else { std::cout << "None" << std::endl; }
    s = onsegment(3.5, c);
    if (s) { (*s).print(); 
    } else { std::cout << "None" << std::endl; }
    s = onsegment(6.5, c);
    if (s) { (*s).print(); 
    } else { std::cout << "None" << std::endl; }
    s = onsegment(0.2, c);
    if (s) { (*s).print(); 
    } else { std::cout << "None" << std::endl; }
    s = onsegment(0.6, c);
    if (s) { (*s).print(); 
    } else { std::cout << "None" << std::endl; }
    s = onsegment(2.5, c);
    if (s) { (*s).print(); 
    } else { std::cout << "None" << std::endl; }
    s = onsegment(4.5, c);
    if (s) { (*s).print(); 
    } else { std::cout << "None" << std::endl; }

    auto as = onantisegment(1.5, c);
    if (as) { (*as).print(); 
    } else { std::cout << "None" << std::endl; }
    as = onantisegment(3.5, c);
    if (as) { (*as).print(); 
    } else { std::cout << "None" << std::endl; }
    as = onantisegment(6.5, c);
    if (as) { (*as).print(); 
    } else { std::cout << "None" << std::endl; }
    as = onantisegment(0.2, c);
    if (as) { (*as).print(); 
    } else { std::cout << "None" << std::endl; }
    as = onantisegment(0.6, c);
    if (as) { (*as).print(); 
    } else { std::cout << "None" << std::endl; }
    as = onantisegment(2.5, c);
    if (as) { (*as).print(); 
    } else { std::cout << "None" << std::endl; }
    as = onantisegment(4.5, c);
    if (as) { (*as).print(); 
    } else { std::cout << "None" << std::endl; }

    
    for (auto idx : nda::range(3) ) {
        auto [i, f] =  segments(c).indices(idx);
        std::cout << i << " " << f << std::endl;
    }

    for (auto idx : nda::range(c.length())) {
        Configuration ctmp = c;
        remove_segment(ctmp, idx);
        for (auto s : segments(ctmp)) { s.print(); }
    }

    for (auto idx : nda::range(c.length())) {
        Configuration ctmp = c;
        remove_antisegment(ctmp, idx);
        for (auto s : antisegments(ctmp)) { s.print(); }
    }

    auto moves = std::vector<MoveFunc>{ NewSegmentInsertionMove(), 
                                        //NewAntiSegmentInsertionMove(), 
                                        //NewSegmentRemoveMove(), 
                                        //NewAntiSegmentRemoveMove() 
    };

    c = Configuration(nda::zeros<double>(0), nda::zeros<double>(0));
    auto S = Solver(c, Delta, e, moves, nt) ;
    S.run();

  return 0;
}
