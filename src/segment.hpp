# pragma once
#include "iostream"
#include <nda/nda.hpp>

#include "configuration.hpp"

struct Segment {
    double t_i;
    double t_f;

    Segment(double i, double f) {
        t_i = i;
        t_f = f;
    }

    double length(double beta) {
        if (t_i < t_f) {
            return t_f - t_i;
        } else {
            return beta - t_i + t_f;
        }
    }

    bool onsegment( double t) {
        if (t_i < t_f ) {
            return (t > t_i) && (t < t_f);
        } else {
            return (t < t_f) || (t > t_i);
        }
    }

    void print() const { std::cout << "Segment(" << t_i << ", "<< t_f << ")" << std::endl; }

};

class SegmentIterator {
    public:
        SegmentIterator(const Configuration& c) : c(c), _state(0) {}

        int length() const { return c.length(); }

        std::pair<int, int> indices(int idx)  {
            if (c.t_f[0] > c.t_i[0]) {
                return std::make_pair(idx, idx);
            } else {
                return std::make_pair(idx, (idx+1) % length());
            }
        }

        Segment getindex(int idx) {
            auto [i_idx, f_idx] = indices(idx);
            return Segment(c.t_i[i_idx], c.t_f[f_idx]);
        }

        SegmentIterator begin() { return *this; }

        SegmentIterator end() {
            SegmentIterator iter = *this;
            iter._state = length();
            return iter;
        }

        SegmentIterator& operator++() {
            _state++;
            return *this;
        }

        SegmentIterator operator++(int){
            SegmentIterator temp = *this;
            ++(*this);
            return temp;
        }

        Segment operator*() {
            int l = length() - 1;
            auto [i_idx, f_idx] = indices(_state);
            if (i_idx  <= l) {
                return Segment(c.t_i[i_idx], c.t_f[f_idx]);
            } else {
                throw std::out_of_range("SegmentIterator out of range");
            }
        }

        bool operator!=(const SegmentIterator& other) const {
            return _state != other._state;
        }

    private:
        const Configuration& c;
        size_t _state;
};

SegmentIterator segments(Configuration& c) { return SegmentIterator(c); }

std::optional<Segment> onsegment(double t, Configuration& c) {
    for (auto s : segments(c) ) {
        if (s.onsegment(t)) { return s; }
    }
    return std::nullopt;
}


void remove_segment(Configuration& c, int segment_idx){
    auto [i_idx, f_idx] = segments(c).indices(segment_idx);
    c.t_i = deleteat(c.t_i, i_idx);
    c.t_f = deleteat(c.t_f, f_idx);
}


bool is_segment_proper(const Configuration& c) {

    if (c.t_i[0] < c.t_f[0]) {
        for (int idx=1; idx < c.length(); idx++) {
            if ( !(c.t_f[idx-1] < c.t_i[idx]) || !(c.t_i[idx] < c.t_f[idx]) ) { return false; }
        }
    } else {
        for (int idx=1; idx < c.length(); idx++) {
            if ( !(c.t_i[idx-1] < c.t_f[idx]) || !(c.t_f[idx] < c.t_i[idx]) ) { return false; }
        }
    }
    return true;
}
