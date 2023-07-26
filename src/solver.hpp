#include <functional>


int randomint(int min, int max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
}

// moves
struct InsertMove {
    double t_i;
    double t_f;
    double l;
    InsertMove(double t_i, double t_f, double l) : t_i(t_i), t_f(t_f), l(l) {}
}

struct RemovalMove {
    double t_i;
    double t_f;
    double l;
    RemovalMove(double t_i, double t_f, double l) : t_i(t_i), t_f(t_f), l(l) {}
}

InsertMove new_segment_insertion_move(const Configuration& c, const Expansion& e) {
    double t_i = e.beta * rand();

    if (c.length() == 0) {
        double t_f = e.beta * nda::rand<>();
        double l = e.beta
    } else {
        auto s = onantisegment(t_i, c);
        double l ( s != std::nulltpr ) ? 0.0 : Segment(t_i, s.t_f).length(e.beta); 
        double t_f = (t_i + l * nda::rand<>() ) % e.beta;
    }

    return InsertMove(t_i, t_f, l);
}

InsertMove new_antisegment_insertion_move(const Configuration& c, const Expansion& e) {
    double t_i = e.beta * nda::rand<>();

    if (c.length() == 0) {
        double t_f = e.beta * nda::rand<>();
        double l = e.beta
    } else {
        auto s = segment(t_i, c);
        double l ( s != std::nulltpr ) ? 0.0 : AntiSegment(t_i, s.t_f).length(e.beta); 
        double t_f = (t_i + l * nda::rand<>() ) % e.beta;
    }

    return InsertMove(t_i, t_f, l);
}

RemovalMove new_segment_remove_move(const Configuration& c, const Expansion& e) {

    if (c.length() > 0) {
        auto idx = randomint(0, c.length());
        auto [i_idx, f_idx] = segments(c).indices(idx);
        auto s = Segment(c.t_i(i_dix), c.t_i( (i_idx+1) % c.length()) );
        auto l = s.length(e.beta);
        return RemovalMove(i_idx, f_idx, l);
    } else {
        return RemovalMove(0.0, 0.0, 0.0);
    }
}

RemovalMove new_antisegment_remove_move(const Configuration& c, const Expansion& e) {

    if (c.length() > 0) {
        auto idx = randomint(0, c.length());
        auto [f_idx, i_idx] = antisegments(c).indices(idx);
        auto s = AntiSegment(c.t_f(f_dix), c.t_f( (f_idx+1) % c.length()) );
        auto l = s.length(e.beta);
        return RemovalMove(i_idx, f_idx, l);
    } else {
        return RemovalMove(0.0, 0.0, 0.0);
    }
}
struct Expansion {
    double beta;
    double h;
    Hybridization Delta;
    std::vector<Move> moves;
    nda::vector<double> move_prop;
    nda::vector<double> move_acc;

    Expansion(double beta, double h, const Hybridization& Delta, std::vector<Move> moves ) 
        : beta(beta), h(h), Delta(Delta), moves(moves),
          move_prop(nda::zeros<double>(move.shape()[0])), 
          move_acc(nda::zeros<double>(moves.shape()[0])) {}
};

struct Determinant {
    private :
        const Configuration c;
        const Expansion e;
        nda::vector<double> t_i;
        nda::vector<double> t_f;
    public:
        nda::matrix<double> mat;
        double value;

        Determinant( Configuation& config, Expansion& expan){
            c = config;
            e = expan;
            t_i = nda::zeros<double>(c.length());
            std::copy(std::begin(c.t_i), std::end(c.t_i), std::begin(t_i));
            t_f = nda::zeros<double>(c.length());
            std::copy(std::begin(c.t_f), std::end(c.t_f), std::begin(t_f));

            if (len(t_f) > 0 && t_f(0) < t_i(0)) {
                t_f = 
            }
        }
};
