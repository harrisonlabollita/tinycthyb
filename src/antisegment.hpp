# pragma once

#include "iostream"
#include <nda/nda.hpp>

#include "configuration.hpp"

struct AntiSegment {
    double t_i;
    double t_f;

    AntiSegment(double t_i, double t_f) : t_i(t_i), t_f(t_f) {}

    double length(const double beta) const {
        if (t_i < t_f) {
            return t_f - t_i;
        } else {
            return beta - t_i + t_f;
        }
    }

    bool onantisegment(const double t) const {
        if (t_i < t_f ) { return (t > t_i) && (t < t_f);
        } else { return (t < t_f) || (t > t_i); }
    }

    void print() const { std::cout << "AntiSegment(" << t_i << ", "<< t_f << ")" << std::endl; }

};

class AntiSegmentIterator {
    public:
        AntiSegmentIterator(const Configuration& configuration) : c(configuration), _state(0) {}

        int length() const { return c.length(); }

        std::pair<int, int> indices(int idx) const {
            if (c.t_f[0] < c.t_i[0]) {
                return std::make_pair(idx, idx);
            } else {
                return std::make_pair(idx, (idx+1) % length());
            }
        }

        AntiSegment getindex(int idx) const {
            auto [f_idx, i_idx] = indices(idx);
            return AntiSegment(c.t_f[f_idx], c.t_i[i_idx]);
        }

        AntiSegmentIterator begin() { return *this; }

        AntiSegmentIterator end() {
            AntiSegmentIterator iter = *this;
            iter._state = length();
            return iter;
        }

        AntiSegmentIterator& operator++() {
            _state++;
            return *this;
        }

        AntiSegmentIterator operator++(int){
            AntiSegmentIterator temp = *this;
            ++(*this);
            return temp;
        }

        AntiSegment operator*() {
            int l = length() - 1;
            auto [f_idx, i_idx] = indices(_state);
            if (f_idx  <= l) { return AntiSegment(c.t_f[f_idx], c.t_i[i_idx]);
            } else { 
                throw std::out_of_range("AntiSegmentIterator out of range");
            }
        }

        bool operator!=(const AntiSegmentIterator& other) const {
            return _state != other._state;
        }

    private:
        const Configuration& c;
        size_t _state;
};

AntiSegmentIterator antisegments(Configuration& c) { return AntiSegmentIterator(c); }

std::optional<AntiSegment> onantisegment(double t, Configuration& c) {
    for (auto s : antisegments(c) ) {
        if (s.onantisegment(t)) { return s; }
    }
    return std::nullopt;
}

void remove_antisegment(Configuration& c, int segment_idx ) {
    auto [f_idx, i_idx] = antisegments(c).indices(segment_idx);
    c.t_i = deleteat(c.t_i, i_idx);
    c.t_f = deleteat(c.t_f, f_idx);
}
