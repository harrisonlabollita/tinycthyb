#pragma once

#include <math.h>
#include <functional>
#include <nda/nda.hpp>
#include <nda/linalg/det_and_inverse.hpp>

#include "segment.hpp"
#include "antisegment.hpp"
#include "util.hpp"

struct Expansion {
    double beta;
    double h;
    Hybridization Delta;
    //std::vector<Move> moves;
    //nda::vector<double> move_prop;
    //nda::vector<double> move_acc;

    Expansion(double beta, double h, Hybridization& Delta) : beta(beta), h(h), Delta(Delta) {}
          //moves(moves),
          //move_prop(nda::zeros<double>(move.size())), 
          //move_acc(nda::zeros<double>(moves.size())) {}
};

struct Determinant {
    public:
        nda::vector<double> t_i;
        nda::vector<double> t_f;
        nda::matrix<double> mat;
        double value;

        Determinant( Configuration& c, Expansion& e){
            t_i = nda::zeros<double>(c.length());
            std::copy(std::begin(c.t_i), std::end(c.t_i), std::begin(t_i));
            t_f = nda::zeros<double>(c.length());
            std::copy(std::begin(c.t_f), std::end(c.t_f), std::begin(t_f));
            if (t_f.size() > 0 && t_f(0) < t_i(0)) {
                t_f = roll(t_f, -1);
            }

            mat = nda::zeros<double>(t_i.size(), t_f.size());
            for (int i = 0; i < t_i.size(); i++) {
                for (int f=0; f < t_f.size(); f++) {
                    mat(f,i) = e.Delta(t_f(f)-t_i(i));
                }
            }

            value = determinant(mat);
        }
};


double trace(Configuration& c, Expansion& e) {
    if (c.length() == 0) { return 2.0; }
    if (!is_segment_proper(c)) { return 0.0; }
    double value = (c.t_f(0) < c.t_f(0)) ? -1.0 : +1.0;

    for (auto s : segments(c)) {
        value *= std::exp(-1*e.h * s.length(e.beta));
    }
    return value;
}

double eval(Configuration& c, Expansion& e) {
    return trace(c, e) * Determinant(c, e).value;
}

// moves
struct InsertMove {
    double t_i;
    double t_f;
    double l;
    InsertMove(double t_i, double t_f, double l) 
        : t_i(t_i), t_f(t_f), l(l) {}
};

struct RemovalMove {
    double t_i;
    double t_f;
    double l;
    RemovalMove(double t_i, double t_f, double l) 
        : t_i(t_i), t_f(t_f), l(l) {}
};

InsertMove new_segment_insertion_move(Configuration& c, Expansion& e) {
    double l;
    double t_f;
    double t_i = e.beta * nda::rand<>();
    if (c.length() == 0) {
        t_f = e.beta * nda::rand<>();
        l   = e.beta;
    } else {
        auto s = onantisegment(t_i, c);
        l = (s) ? Segment(t_i, (*s).t_f).length(e.beta) : 0.0;
        t_f = fmod( (t_i + l * nda::rand<>() ), e.beta);
    }
    return InsertMove(t_i, t_f, l);
}


InsertMove new_antisegment_insertion_move(Configuration& c, Expansion& e) {
    double l;
    double t_f;
    double t_i = e.beta * nda::rand<>();
    if (c.length() == 0) {
        t_f = e.beta * nda::rand<>();
        l = e.beta;
    } else {
        auto s = onsegment(t_i, c);
        l = (s) ? AntiSegment(t_i, (*s).t_f).length(e.beta) : 0.0;
        t_f = fmod((t_i + l * nda::rand<>() ), e.beta);
    }

    return InsertMove(t_i, t_f, l);
}

RemovalMove new_segment_remove_move(Configuration& c, Expansion& e) {

    if (c.length() > 0) {
        auto idx = randomint(0, c.length());
        auto [i_idx, f_idx] = segments(c).indices(idx);
        auto s = Segment(c.t_i(i_idx), c.t_i( (i_idx+1) % c.length()) );
        auto l = s.length(e.beta);
        return RemovalMove(i_idx, f_idx, l);
    } else {
        return RemovalMove(0.0, 0.0, 0.0);
    }
}

RemovalMove new_antisegment_remove_move(Configuration& c, Expansion& e) {

    if (c.length() > 0) {
        auto idx = randomint(0, c.length());
        auto [f_idx, i_idx] = antisegments(c).indices(idx);
        auto s = AntiSegment(c.t_f(f_idx), c.t_f( (f_idx+1) % c.length()) );
        auto l = s.length(e.beta);
        return RemovalMove(i_idx, f_idx, l);
    } else {
        return RemovalMove(0.0, 0.0, 0.0);
    }
}
