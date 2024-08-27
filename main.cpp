/*
    Args:
    M: large particle mass
    m: small particle mass
    N_width: width (in particle number) of distribution
    N_length: length (in particle number) of distribution
    d: distance between particles
    v0: initial speed of large particle
    D: initial distance of large particle from distribution (negative is inside)
    T: integration time
    N: integration steps
    file: where to store results
*/

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cmath>

const double G = 6.67e-11;
using vec = std::vector<double>;


struct solver {
    double _M;
    double _m;
    double _step;
    unsigned _N_steps;
    std::vector<vec>& _positions;
    std::vector<vec>& _velocities;
    vec accels;
    
    // positions has size 3*(N_particles + 1)
    // positions = [x1, y1, z1, x2, y2, z2, ...]
    vec calc_accelerations(const vec& positions) {
        vec accels(positions.size(), 0);
        std::size_t N_parts = positions.size()/3;
        for (std::size_t i = 1; i < N_parts; i++) {
            std::array<double, 3> pos_rel {
                positions[3*i] - positions[0],
                positions[3*i+1] - positions[1],
                positions[3*i+2] - positions[2]
            };
            double inv_d3 = G / std::pow(std::sqrt(pos_rel[0]*pos_rel[0]+pos_rel[1]*pos_rel[1]+pos_rel[2]*pos_rel[2]), 3);
            for (int j : {0, 1, 2}) {
                accels[j] += _m * inv_d3 * pos_rel[j];
                accels[3*i+j] -= _M * inv_d3 * pos_rel[j];
            }
        }
        return accels;
    }

    void multiply(vec& v, double a) {
        for (auto& it : v) {
            it *= a;
        }
    }

    void add_to_first(vec& v1, const vec& v2) {
        for (unsigned i = 0; i < v1.size(); i++) {
            v1[i] += v2[i];
        }
    }

    void solve() {
        for (unsigned i = 0; i < _positions.size()-1; i++) {
            static unsigned print_interval = (unsigned) _positions.size()/100;
            if (i % print_interval == 0) {
                std::cout << "\rt = " << i*_step << std::flush;
            }
            vec& pos = _positions[i];
            vec& vel = _velocities[i];

            vec k1x = vel;
            multiply(k1x, _step/2);
            vec k1v = calc_accelerations(pos); // Should use move semantics
            multiply(k1v, _step/2);
            
            vec k2x = vel;
            add_to_first(k2x, k1v);
            multiply(k2x, _step/2);
            vec k2v = pos;
            add_to_first(k2v, k1x);
            k2v = calc_accelerations(k2v);
            multiply(k2v, _step/2);

            vec k3x = vel;
            add_to_first(k3x, k2v);
            multiply(k3x, _step);
            vec k3v = pos;
            add_to_first(k3v, k2x);
            k3v = calc_accelerations(k3v);
            multiply(k3v, _step);

            vec k4x = vel;
            add_to_first(k4x, k3v);
            multiply(k4x, _step);
            vec k4v = pos;
            add_to_first(k4v, k3x);
            k4v = calc_accelerations(k4v);
            multiply(k4v, _step);

            multiply(k1x, 1.0/3);
            multiply(k1v, 1.0/3);
            multiply(k2x, 2.0/3);
            multiply(k2v, 2.0/3);
            multiply(k3x, 2.0/3);
            multiply(k3v, 2.0/3);
            multiply(k4x, 1.0/6);
            multiply(k4v, 1.0/6);
            add_to_first(k3x, k4x);
            add_to_first(k3v, k4v);
            add_to_first(k2x, k3x);
            add_to_first(k2v, k3v);
            add_to_first(k1x, k2x);
            add_to_first(k1v, k2v);

            _positions[i+1] = pos;
            _velocities[i+1] = vel;
            add_to_first(_positions[i+1], k1x);
            add_to_first(_velocities[i+1], k1v);
        }
        std::cout << "\n";
    }
};

int main(int argc, char* argv[]) {
    if (argc != 11) {
        std::cout << "Received " << argc-1 << " arguments." << "\n";
        return 1;
    }

    double M = std::stod(argv[1]);
    double m = std::stod(argv[2]);
    unsigned N_width = std::stoi(argv[3]);
    unsigned N_length = std::stoi(argv[4]);
    double distance = std::stod(argv[5]);
    double v0 = std::stod(argv[6]);
    double d0 = std::stod(argv[7]);
    double T = std::stod(argv[8]);
    unsigned N = std::stoi(argv[9]);
    std::string filename = argv[10];

    unsigned N_particles = 2*N_width * 2*N_width * N_length;
    std::vector<vec> positions(N+1, vec(3*(N_particles+1), 0));
    std::vector<vec> velocities(N+1, vec(3*(N_particles+1), 0));

    auto& pos_init = positions[0];
    pos_init[0] = -d0;
    double x_min = -distance/2 - N_width * distance;
    double y_min = x_min;
    for (unsigned k = 0; k < N_length; k++) {
        for (unsigned i = 0; i < 2*N_width; i++) {
            for (unsigned j = 0; j < 2*N_width; j++) {
                pos_init[3+i+j+k] = x_min + i*distance;
                pos_init[3+i+j+k+1] = y_min + j*distance;
                pos_init[3+i+j+k+2] = k*distance;
            }
        }
    }
    velocities[0][0] = v0;

    solver Solver {M, m, T/N, N, positions, velocities, vec(3*(N_particles+1), 0)};
    Solver.solve();

    return 0;
}