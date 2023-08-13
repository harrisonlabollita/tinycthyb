#pragma once

#include <functional>
#include <limits>
#include <math.h>
#include <nda/linalg/det_and_inverse.hpp>
#include <nda/nda.hpp>
#include <variant>
#include <vector>

#include "configuration.hpp"
#include "segment.hpp"
#include "antisegment.hpp"
#include "hybridization.hpp"
#include "util.hpp"

namespace tinycthyb {

struct Expansion {
  double beta;
  double h;
  Hybridization &Delta;
  Expansion(double beta, double h, Hybridization &Delta)
      : beta(beta), h(h), Delta(Delta) {}
};

struct Determinant {
public:
  nda::vector<double> t_i;
  nda::vector<double> t_f;
  nda::matrix<double> mat;
  double value;

  Determinant(Configuration &c, Expansion &e) {
    t_f = c.t_f;
    t_i = c.t_i;

    if (t_f.size() > 0 && t_f(0) < t_i(0)) {
      t_f = roll(t_f, -1);
    }

    mat = nda::zeros<double>(t_i.size(), t_f.size());
    for (int i = 0; i < t_i.size(); i++) {
      for (int j = 0; j < t_f.size(); j++) {
        mat(j, i) = e.Delta(t_f(j) - t_i(i));
      }
    }
    value = determinant(mat);
  }
};

double trace(Configuration &c, Expansion &e) {
  if (c.length() == 0) {
    return 2.0;
  }
  if (!is_segment_proper(c)) {
    return 0.0;
  }
  double value = (c.t_f(0) < c.t_i(0)) ? -1.0 : +1.0;

  for (auto s : segments(c)) {
    value *= std::exp(-e.h * s.length(e.beta));
  }
  return value;
}

double eval(Configuration &c, Expansion &e) {
  return trace(c, e) * Determinant(c, e).value;
}

// moves
InsertMove NewSegmentInsertionMove(Configuration &c, Expansion &e) {
  double l;
  double t_f;
  double t_i = e.beta * nda::rand<>();
  if (c.length() == 0) {
    t_f = e.beta * nda::rand<>();
    l = e.beta;
  } else {
    auto s = onantisegment(t_i, c);
    if (s) {
      l = Segment(t_i, (*s).t_f).length(e.beta);
    } else {
      l = 0.0;
    }
    t_f = fmod((t_i + l * nda::rand<>()), e.beta);
  }
  return InsertMove(t_i, t_f, l);
}

InsertMove NewAntiSegmentInsertionMove(Configuration &c, Expansion &e) {
  double l;
  double t_f;
  double t_i = e.beta * nda::rand<>();
  if (c.length() == 0) {
    t_f = e.beta * nda::rand<>();
    l = e.beta;
  } else {
    auto s = onsegment(t_i, c);
    if (s) {
      l = AntiSegment(t_i, (*s).t_f).length(e.beta);
    } else {
      l = 0.0;
    }
    t_f = fmod((t_i + l * nda::rand<>()), e.beta);
  }
  return InsertMove(t_f, t_i, l);
}

RemovalMove NewSegmentRemoveMove(Configuration &c, Expansion &e) {
  if (c.length() > 0) {
    auto idx = randomint(0, c.length() - 1);
    auto [i_idx, f_idx] = segments(c).indices(idx);
    auto s = Segment(c.t_i(i_idx), c.t_i((i_idx + 1) % c.length()));
    auto l = s.length(e.beta);
    return RemovalMove(i_idx, f_idx, l);
  } else {
    return RemovalMove(0, 0, 0.0);
  }
}

RemovalMove NewAntiSegmentRemoveMove(Configuration &c, Expansion &e) {
  if (c.length() > 0) {
    auto idx = randomint(0, c.length() - 1);
    auto [f_idx, i_idx] = antisegments(c).indices(idx);
    auto s = AntiSegment(c.t_f(f_idx), c.t_f((f_idx + 1) % c.length()));
    auto l = s.length(e.beta);
    return RemovalMove(i_idx, f_idx, l);
  } else {
    return RemovalMove(0, 0, 0.0);
  }
}

using Moves = std::variant<InsertMove, RemovalMove>;
using MoveFunc = std::function<Moves(Configuration &, Expansion &)>;

class Solver {
private:
  Hybridization &Delta;
  Expansion &e;

public:
  std::vector<MoveFunc> moves;
  GreensFunction g;
  int nt;
  nda::vector<double> move_prop;
  nda::vector<double> move_acc;

  Solver(Hybridization &Delta, Expansion &e, std::vector<MoveFunc> moves,
         int nt)
      : Delta(Delta), e(e), moves(moves), nt(nt), g(e.beta, nt) {
    move_prop = nda::zeros<double>(moves.size());
    move_acc = nda::zeros<double>(moves.size());
  }

  void sample_greens_function(Configuration &c) {
    auto d = Determinant(c, e);
    auto M = inverse(d.mat);
    auto w = trace(c, e) * d.value;
    g.sign += sign(w);
    for (auto i = 0; i < c.length(); i++) {
      for (auto j = 0; j < c.length(); j++) {
        g.accumulate(d.t_f(i) - d.t_i(j), M(j, i));
      }
    }
  }

  double propose(Configuration &c, InsertMove &move) {
    auto cnew = c + move;
    double R = move.l * e.beta / cnew.length() *
               (std::abs(eval(cnew, e) / eval(c, e)));
    return R;
  }

  double propose(Configuration &c, RemovalMove &move) {
    auto cnew = c + move;
    if (move.l == 0) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    double R =
        c.length() / e.beta / move.l * (std::abs(eval(cnew, e) / eval(c, e)));
    return R;
  }

  Configuration finalize(Configuration &c, InsertMove &move) {
    auto cnew = c + move;
    return cnew;
  }

  Configuration finalize(Configuration &c, RemovalMove &move) {
    auto cnew = c + move;
    return cnew;
  }

  Configuration metropolis_hastings_update(Configuration &c) {
    auto move_idx = randomint(0, moves.size() - 1);
    auto m = moves[move_idx](c, e);
    double R = 0.0;
    if (std::holds_alternative<InsertMove>(m)) {
      InsertMove move = std::get<InsertMove>(m);
      move_prop(move_idx) += 1;
      R = propose(c, move);
      if (R > nda::rand<>()) {
        c = finalize(c, move);
        move_acc(move_idx) += 1;
      }
    } else if (std::holds_alternative<RemovalMove>(m)) {
      RemovalMove move = std::get<RemovalMove>(m);
      move_prop(move_idx) += 1;
      R = propose(c, move);
      if (R > nda::rand<>()) {
        c = finalize(c, move);
        move_acc(move_idx) += 1;
      }
    }
    return c;
  }

  void solve(Configuration c, int epoch_steps = 10, int warmup_epochs = 1000,
             long sampling_epochs = 100000) {

    std::cout << "Starting CT-HYB QMC" << std::endl;

    std::cout << "Warmup epochs " << warmup_epochs << " with " << epoch_steps
              << " steps." << std::endl;

    for (auto epoch = 0; epoch < warmup_epochs; epoch++) {
      for (auto step = 0; step < epoch_steps; step++) {
        c = metropolis_hastings_update(c);
      }
    }

    std::cout << "Sampling epochs " << sampling_epochs << " with "
              << epoch_steps << " steps." << std::endl;

    for (auto epoch = 0; epoch < sampling_epochs; epoch++) {
      for (auto step = 0; step < epoch_steps; step++) {
        c = metropolis_hastings_update(c);
      }
      sample_greens_function(c);
    }
  }
};

} // namespace tinycthyb
