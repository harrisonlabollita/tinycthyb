#pragma once

#include <variant>
#include <vector>
#include <math.h>
#include <limits>
#include <functional>
#include <nda/nda.hpp>
#include <nda/linalg/det_and_inverse.hpp>

#include "configuration.hpp"
#include "segment.hpp"
#include "antisegment.hpp"
#include "util.hpp"

struct Expansion {
    double beta;
    double h;
    Hybridization Delta;
    Expansion(double beta, double h, Hybridization& Delta) : beta(beta), h(h), Delta(Delta) {}
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
struct NewSegmentInsertionMove {
    InsertMove operator()(Configuration& c, Expansion& e) {
        std::cout << "called 1" << std::endl;
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
        std::cout << t_i << " "<< t_f << " "<< l << std::endl;
        return InsertMove(t_i, t_f, l);
    }
};



struct NewAntiSegmentInsertionMove {
    InsertMove operator()(Configuration& c, Expansion& e) {
        std::cout << "called 2" << std::endl;
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

        std::cout << t_i << " "<< t_f << " "<< l << std::endl;
        return InsertMove(t_i, t_f, l);
    }
};

struct NewSegmentRemoveMove {
    RemovalMove operator()(Configuration& c, Expansion& e) {
        std::cout << "called 3" << std::endl;
        if (c.length() > 0) {
            auto idx = randomint(0, c.length());
            auto [i_idx, f_idx] = segments(c).indices(idx);
            auto s = Segment(c.t_i(i_idx), c.t_i( (i_idx+1) % c.length()) );
            auto l = s.length(e.beta);
            return RemovalMove(i_idx, f_idx, l);
            std::cout << i_idx << " " << f_idx << " " << l << std::endl;
        } else {
            return RemovalMove(0, 0, 0.0);
        }
    }
};


struct NewAntiSegmentRemoveMove {
    RemovalMove operator()(Configuration& c, Expansion& e) {
        std::cout << "called 4" << std::endl;
        if (c.length() > 0) {
            auto idx = randomint(0, c.length());
            std::cout << idx << std::endl;
            auto [f_idx, i_idx] = antisegments(c).indices(idx);
            std::cout << f_idx<< std::endl;
            std::cout << i_idx<< std::endl;
            auto s = AntiSegment(c.t_f(f_idx), c.t_f( (f_idx+1) % c.length()) );
            s.print();
            auto l = s.length(e.beta);
            std::cout << l << std::endl;
            return RemovalMove(i_idx, f_idx, l);
        } else {
            return RemovalMove(0, 0, 0.0);
        }
    }
};

using Moves    = std::variant<InsertMove, RemovalMove>;
using MoveFunc = std::function<Moves(Configuration&, Expansion&)>;


class Solver {
    private:
        Configuration c;
        Hybridization Delta;
        Expansion     e;
        std::vector<MoveFunc> moves;
    public:
        GreensFunction g;
        int nt;
        nda::vector<double> move_prop;
        nda::vector<double> move_acc;

        Solver(Configuration& c, Hybridization& Delta, Expansion& e, std::vector<MoveFunc> moves, int nt ) 
            : c(c), Delta(Delta), e(e), moves(moves), nt(nt), g(e.beta, nt) {
            move_prop = nda::zeros<double>(moves.size());
            move_acc  = nda::zeros<double>(moves.size());
        }

        void sample_greens_function() {
            auto  d = Determinant(c, e);
            auto M = inverse(d.mat);

            auto w = trace(c, e) * d.value;

            g.sign += sign(w);

            for (auto i : nda::range(nt)) {
                for (auto j : nda::range(nt)) { g.accumulate(d.t_f(i)-d.t_i(j), M(j,i));}
            }
        }

        double propose(InsertMove& move) {
            auto cnew = c + move;
            double R = move.l * e.beta / cnew.length() * ( std::abs(eval(cnew, e) / eval(c, e) ) );
            return R;
        }

        double propose(RemovalMove& move) {
            if (move.l == 0) { return std::numeric_limits<double>::quiet_NaN(); }
            auto cnew = c + move;
            double R = c.length() / e.beta / move.l * (std::abs(eval(cnew, e) / eval(c, e) ));
            return R;
        }

        void finalize(InsertMove& move) {
            c = c + move;
        }

        void finalize(RemovalMove& move) {
            c = c + move;
        }

        void metropolis_hastings_update() {
            c.print();
            auto move_idx = randomint(0, moves.size());
            std::cout << move_idx << std::endl; 
            auto m = moves[move_idx](c, e);
            if (std::holds_alternative<InsertMove>(m))  {
                InsertMove move = std::get<InsertMove>(m);
                move_prop[move_idx] += 1;
                auto R =  propose(move);
                std::cout << "R = " << R << std::endl;
                if (R > nda::rand<>()) {
                    finalize(move);
                    move_acc[move_idx] += 1;
                }
            } else if (std::holds_alternative<RemovalMove>(m) ) {
                RemovalMove move = std::get<RemovalMove>(m);
                move_prop[move_idx] += 1;
                auto R =  propose(move);
                std::cout << "R = " << R << std::endl;
                if (R > nda::rand<>()) {
                    finalize(move);
                    move_acc[move_idx] += 1;
                }
            }
        }

        void run(int epoch_steps=10, int warmup_epochs=1, long sampling_epochs = 1) {

            std::cout << "Starting CT-HYB QMC" << std::endl;

            std::cout << "Warmup epochs " << warmup_epochs << " with " << epoch_steps << " steps." << std::endl;

            for (auto epoch : nda::range(warmup_epochs)) {
                for (auto step : nda::range(epoch_steps)) { metropolis_hastings_update(); }
            }

            for (auto epoch : nda::range(sampling_epochs)) {
                for (auto step : nda::range(epoch_steps)) { metropolis_hastings_update(); }
                sample_greens_function();
            }
        }
};
