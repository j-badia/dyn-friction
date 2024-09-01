/*
    Args:
    M: large particle mass
    m: small particle mass
    N_width: width (in particle number) of distribution
    N_length: length (in particle number) of distribution
    d: distance between particles
    v0: initial speed of large particle
    D: initial distance of large particle from distribution (positive is inside)
    T: integration time
    N: integration steps
    file: where to store results
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <chrono>

const std::string default_params_file {"params"};
const double G = 4.3e-3; // In pc (km/s)^2 / M_sun
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
            double d2 = pos_rel[0]*pos_rel[0]+pos_rel[1]*pos_rel[1]+pos_rel[2]*pos_rel[2];
            double inv_d3 = G / std::sqrt(d2*d2*d2);
            //double inv_d3 = G / std::pow(std::sqrt(pos_rel[0]*pos_rel[0]+pos_rel[1]*pos_rel[1]+pos_rel[2]*pos_rel[2]), 3);
            //double inv_d3 = G * std::pow(pos_rel[0]*pos_rel[0]+pos_rel[1]*pos_rel[1]+pos_rel[2]*pos_rel[2], -1.5);
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
        unsigned N = v1.size();
        for (unsigned i = 0; i < N; i++) {
            v1[i] += v2[i];
        }
    }

    void add_to_first_then_mult(vec& v1, const vec& v2, const double a) {
        unsigned N = v1.size();
        for (unsigned i = 0; i < N; i++) {
            v1[i] += v2[i];
            v1[i] *= a;
        }
    }

    void mult_then_add_to_first(vec& v1, const vec& v2, const double a) {
        unsigned N = v1.size();
        for (unsigned i = 0; i < N; i++) {
            v1[i] += a * v2[i];
        }
    }

    void solve() {
        for (unsigned i = 0; i < _positions.size()-1; i++) {
            static unsigned print_interval = (unsigned) _positions.size()/100;
            if (i % print_interval == 0) {
                std::cout << "\rt = " << i*_step << std::flush;
            }
            const vec& pos = _positions[i];
            const vec& vel = _velocities[i];

            vec k1x = vel;
            multiply(k1x, _step/2);
            vec k1v = calc_accelerations(pos); // Should use move semantics
            multiply(k1v, _step/2);
            
            vec k2x = vel;
            add_to_first_then_mult(k2x, k1v, _step/2);
            vec k2v = pos;
            add_to_first(k2v, k1x);
            k2v = calc_accelerations(k2v);
            multiply(k2v, _step/2);

            vec k3x = vel;
            add_to_first_then_mult(k2x, k2v, _step);
            vec k3v = pos;
            add_to_first(k3v, k2x);
            k3v = calc_accelerations(k3v);
            multiply(k3v, _step);

            vec k4x = vel;
            add_to_first_then_mult(k4x, k3v, _step);
            vec k4v = pos;
            add_to_first(k4v, k3x);
            k4v = calc_accelerations(k4v);
            multiply(k4v, _step);

            _positions[i+1] = pos;
            _velocities[i+1] = vel;
            for (unsigned j = 0; j < pos.size(); j++) {
                _positions[i+1][j] += k1x[j]/3 + 2*k2x[j]/3 + k3x[j]/3 + k4x[j]/6;
                _velocities[i+1][j] += k1v[j]/3 + 2*k2v[j]/3 + k3v[j]/3 + k4v[j]/6;
            }
        }
        std::cout << "\n";
    }
};

int main(int argc, char* argv[]) {
    double M, m, distance, v0, d0, T;
    unsigned N_width, N_length, N;
    std::string filename;
    if (argc == 1 || argc == 2) {
        std::string params_file;
        if (argc == 1) {
            params_file = default_params_file;
        } else if (argc == 2) {
            params_file = argv[1];
        }
        std::ifstream file;
        file.open(params_file);
        std::string name, value;
        while (file.peek() != EOF) {
            std::getline(file, name, '=');
            std::getline(file, value);
            if (name == "M") {
                M = std::stod(value);
            } else if (name == "m") {
                m = std::stod(value);
            } else if (name == "Nw") {
                N_width = std::stoi(value);
            } else if (name == "Nl") {
                N_length = std::stoi(value);
            } else if (name == "D") {
                distance = std::stod(value);
            } else if (name == "v0") {
                v0 = std::stod(value);
            } else if (name == "d0") {
                d0 = std::stod(value);
            } else if (name == "T") {
                T = std::stod(value);
            } else if (name == "N") {
                N = std::stoi(value);
            } else if (name == "filename") {
                filename = value;
            }
        }
    } else if (argc == 11) {
        M = std::stod(argv[1]);
        m = std::stod(argv[2]);
        N_width = std::stoi(argv[3]);
        N_length = std::stoi(argv[4]);
        distance = std::stod(argv[5]);
        v0 = std::stod(argv[6]);
        d0 = std::stod(argv[7]);
        T = std::stod(argv[8]);
        N = std::stoi(argv[9]);
        filename = argv[10];
    } else {
        std::cout << "Received " << argc-1 << " arguments." << "\n";
        return 1;
    }    

    unsigned N_particles = 2*N_width * 2*N_width * N_length;
    std::vector<vec> positions(N+1, vec(3*(N_particles+1), 0));
    std::vector<vec> velocities(N+1, vec(3*(N_particles+1), 0));

    auto& pos_init = positions[0];
    pos_init[2] = d0;
    double x_min = -distance/2 - (N_width-1) * distance;
    double y_min = x_min;
    for (unsigned k = 0; k < N_length; k++) {
        for (unsigned i = 0; i < 2*N_width; i++) {
            for (unsigned j = 0; j < 2*N_width; j++) {
                pos_init[3 + 3*(i + 2*N_width*j + 2*N_width*2*N_width*k)] = x_min + i*distance;
                pos_init[3 + 3*(i + 2*N_width*j + 2*N_width*2*N_width*k) + 1] = y_min + j*distance;
                pos_init[3 + 3*(i + 2*N_width*j + 2*N_width*2*N_width*k) + 2] = k*distance;
            }
        }
    }
    velocities[0][2] = v0;

    solver Solver {M, m, T/N, N, positions, velocities, vec(3*(N_particles+1), 0)};

    auto t0 = std::chrono::steady_clock::now();
    Solver.solve();
    auto t1 = std::chrono::steady_clock::now();
    std::cout << "Time: " << std::chrono::duration<double>{t1-t0}.count() << "s" << std::endl;

    std::ofstream file;
    file.open(filename);
    for (unsigned i = 0; i < positions.size(); i++) {
        file << i*T/N;
        for (int j : {0, 1, 2}) {
            file << " " << positions[i][j];
        }
        for (int j : {0, 1, 2}) {
            file << " " << velocities[i][j];
        }
        file << "\n";
    }
    file.close();

    return 0;
}