/*
    Args:
    M: large particle mass
    m: small particle mass
    N_width: width (in particle number) of distribution
    N_length: length (in particle number) of distribution
    D: distance between particles
    v0: initial speed of large particle
    d0: initial distance of large particle from distribution (positive is inside)
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
using vec3 = std::vector<std::array<double, 3>>;


class Solver {
    /*
    Arrays needed:
    All positions and velocities for the current step: sizes 3*(N_small_particles+1) for both. Will be managed by the class so N_small_particles needs to be a parameter
    All positions and velocities for the next step
    Positions and velocities for the big particle: sizes (N_steps x 3) for both. Caller must allocate, passed by reference.
    */
    private:

    double M;
    double m;
    double step;
    unsigned N_small;
    unsigned N_steps;
    
    vec curr_positions;
    vec curr_velocities;
    vec next_positions;
    vec next_velocities;
    vec accels;

    vec3& big_positions;
    vec3& big_velocities;
    
    // positions has size 3*(N_particles + 1)
    // positions = [x1, y1, z1, x2, y2, z2, ...]
    void calc_accelerations(const vec& positions) {
        accels.assign(accels.size(), 0);
        for (unsigned i = 1; i < N_small+1; i++) {
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
                accels[j] += m * inv_d3 * pos_rel[j];
                accels[3*i+j] -= M * inv_d3 * pos_rel[j];
            }
        }
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

    public:

    Solver(double _M, double _m, double _step, unsigned _N_small, unsigned _N_steps, const vec& init_pos, const vec& init_vel, vec3& _big_pos, vec3& _big_vel)
        : M{_M}, m{_m}, step{_step}, N_small{_N_small}, N_steps{_N_steps}, big_positions{_big_pos}, big_velocities{_big_vel}
    {
        accels.assign(3*(N_small+1), 0);
        curr_positions = init_pos;
        curr_velocities = init_vel;
        next_positions.reserve(3*(N_small+1));
        next_velocities.reserve(3*(N_small+1));

        for (int i : {0, 1, 2}) {
            big_positions[0][i] = init_pos[i];
            big_velocities[0][i] = init_vel[i];
        }
    }

    void solve() {
        unsigned print_interval = (unsigned) N_steps/100;
        for (unsigned i = 1; i < N_steps; i++) {
            if (i % print_interval == 0) {
                std::cout << "\rt = " << i*step << std::flush;
            }

            vec k1x = curr_velocities;
            multiply(k1x, step/2);
            calc_accelerations(curr_positions);
            vec k1v = accels;
            multiply(k1v, step/2);
            
            vec k2x = curr_velocities;
            add_to_first_then_mult(k2x, k1v, step/2);
            vec k2v = curr_positions;
            add_to_first(k2v, k1x);
            calc_accelerations(k2v);
            k2v = accels;
            multiply(k2v, step/2);

            vec k3x = curr_velocities;
            add_to_first_then_mult(k2x, k2v, step);
            vec k3v = curr_positions;
            add_to_first(k3v, k2x);
            calc_accelerations(k3v);
            k3v = accels;
            multiply(k3v, step);

            vec k4x = curr_velocities;
            add_to_first_then_mult(k4x, k3v, step);
            vec k4v = curr_positions;
            add_to_first(k4v, k3x);
            calc_accelerations(k4v);
            k4v = accels;
            multiply(k4v, step);

            next_positions = curr_positions;
            next_velocities = curr_velocities;
            for (unsigned j = 0; j < next_positions.size(); j++) {
                next_positions[j] += k1x[j]/3 + 2*k2x[j]/3 + k3x[j]/3 + k4x[j]/6;
                next_velocities[j] += k1v[j]/3 + 2*k2v[j]/3 + k3v[j]/3 + k4v[j]/6;
            }

            for (int j : {0, 1, 2}) {
                big_positions[i][j] = next_positions[j];
                big_velocities[i][j] = next_velocities[j];
            }

            std::swap(curr_positions, next_positions);
            std::swap(curr_velocities, next_velocities);
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

    unsigned N_small = 2*N_width * 2*N_width * N_length;

    vec pos_init;
    pos_init.assign(3*(N_small+1), 0);
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

    vec vel_init(3*(N_small+1), 0);
    vel_init[2] = v0;

    vec3 big_positions;
    vec3 big_velocities;
    big_positions.reserve(N);
    big_velocities.reserve(N);

    Solver solver {M, m, T/N, N_small, N, pos_init, vel_init, big_positions, big_velocities};

    auto t0 = std::chrono::steady_clock::now();
    solver.solve();
    auto t1 = std::chrono::steady_clock::now();
    std::cout << "Time: " << std::chrono::duration<double>{t1-t0}.count() << "s" << std::endl;

    std::ofstream file;
    file.open(filename);
    for (unsigned i = 0; i < N; i++) {
        file << i*T/N;
        for (int j : {0, 1, 2}) {
            file << " " << big_positions[i][j];
        }
        for (int j : {0, 1, 2}) {
            file << " " << big_velocities[i][j];
        }
        file << "\n";
    }
    file.close();

    return 0;
}