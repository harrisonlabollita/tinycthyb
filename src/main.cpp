// Segment picture hybridization expansion
// Author: H. LaBollita (2023)
//
// C++ (TRIQS/nda) version of https://github.com/HugoStrand/cthyb.jl/blob/main/cthyb.jl
//
//

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

    int length(Configuration& config) { return config.t_i.shape()[0]; }

    void print(Configuration& config) {
        std::cout << "Configuration(" << config.t_i << ", " << config.t_f << ")" << std::endl;
    }

};

int length(Configuration c) { return c.t_i.shape()[0] };

struct Segment {
    double t_i;
    double t_f;
    Segment(double i, double f) {
        t_i = i;
        t_f = f;
    }

    double length(Segment& seg, double beta) {
        if (seg.t_i < seg.t_f) {
            return seg.t_f - seg.t_i
        } else {
            return beta - seg.t_i + seg.t_f
        }
    }

    bool onsegment(Segment& seg, double t) {
        if (seg.t_i < seg.t_f ) {
            return (t > seg.t_i) && (t < seg.t_f)
        } else {
            return (t < seg.t_f) || (t > seg.t_i)
        }
    }
};

// TODO : implement SegmentIterator class
struct SegmentIterator {
    Configuration c;
    SegmentIterator(Configuration config) { c = config };
};


SegmentIterator segments(Configuration c) { return SegmentIterator(c) };


// Anti-segment utilities

struct AntiSegment {
    double t_i;
    double t_f;
    AntiSegment(double i, double f) {
        t_i = i;
        t_f = f;
    }

    double length(AntiSegment& s, double beta) {
        return s.t_i < s.t_f ? s.t_f - s.t_i : beta - s.t_i + s.t_f;
    }

};


struct AntiSegmentIterator {
    Configuration c;
    AntiSegmentIterator(Configuration config){
        c = config;
    }
};

AntiSegmentIterator antisegments(Configuration c) { return AntiSegmentIterator(c); }



void remove_segment(Configuration &c, int segment_idx) {
    auto i_idx, f_idx = indices(segments(c), segment_idx);
    c.t_i = deleteat(c.t_i, i_idx);
    c.t_f = deleteat(c.t_f, f_idx);
};

void remove_antisegment(Configuration &c, int segment_idx) {
    auto f_idx, i_idx = indices(antisegments(c), segment_idx);
    c.t_i = deleteat(c.t_i, i_idx);
    c.t_f = deleteat(c.t_f, f_idx);
};


int main(){
    auto c1 = Configuration(nda::vector<double>{2.0, 1.0},
                           nda::vector<double>{3.0, 4.0}
                           );

    auto c2 = Configuration(nda::vector<double>{},
                           nda::vector<double>{}
                           );

    std::cout << "c.t_i =" << c1.t_i << std::endl;
    std::cout << "c.t_f =" << c1.t_f << std::endl;
    std::cout << "len(c) = " << length(c1) << std::endl;
    std::cout << "c.t_i =" << c2.t_i << std::endl;
    std::cout << "c.t_f =" << c2.t_f << std::endl;
    std::cout << "len(c) = " << length(c2) << std::endl;
  return 0;
}
