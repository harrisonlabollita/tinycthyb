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

        //static_assert(t_i.shape()[0] == t_f.shape()[0]);
        
        if (i.shape()[0] == 0) {
            t_i = nda::vector<double>{};
            t_f = nda::vector<double>{};
        }
        std::sort(i.begin(), i.end());
        std::sort(f.begin(), f.end());
        t_i = i;
        t_f = f;
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
};

double length(Segment s, double beta){return s.t_i < s.t_f ? s.t_f - s.t_i : beta - s.t_i + s.t_f };

//class SegmentIterator {
 //   private: 
  //      nda::vector<double>& 
//}

struct SegmentIterator {
    Configuration c;
    SegmentIterator(Configuration config) { c = config };
};

int length(SegmentIterator s) { return length(s.c) };

SegmentIterator segments(Configuration c) { return SegmentIterator(c) };

// TODO: implement mod
std::tuple<int, int> indices(SegmentIterator, int idx) {
    return s.c.t_f[0] > s.c.t_i[0] ? (idx, idx) : (idx, mod(idx, range(length(s.c));
}

Segment getindex(SegmentIterator s, int idx){
    auto i_idx, f_idx = indices(s, idx);
    return Segment(s.c.t_i[i_idx], s.c.t_f[f_idx]);
}

// Anti-segment utilities

struct AntiSegment {
    double t_i;
    double t_f;
    AntiSegment(double i, double f) {
        t_i = i;
        t_f = f;
    }
};

double length(AntiSegment s, double beta) {
    return s.t_i < s.t_f ? s.t_f - s.t_i : beta - s.t_i + s.t_f;
}

struct AntiSegmentIterator {
    Configuration c;
    AntiSegmentIterator(Configuration config){
        c = config;
    }
};

AntiSegmentIterator antisegments(Configuration c) { return AntiSegmentIterator(c); }


bool onsegment(double t, Segment s) {
    if (s.t_i < s.t_f) {
        return (t > s.t_i) && (t < s.t_f);
    } else {
        return (t < s.t_f) || (t > s.t_i);
    }
}

bool onantisegment(double t, AntiSegment s) {
    if(s.t_i < s.t_f){
        return (t > s.t_i) && (t < s.t_f);
    } else {
        return (t < s.t_f) || (t > s.t_i);
    }
}


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
