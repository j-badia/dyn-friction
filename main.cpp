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
#include <chrono>

#include "solver.hpp"

const std::string default_params_file {"params"};

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
            } else if (name == "file") {
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

    vec3 big_positions(N, {0, 0, 0});
    vec3 big_velocities(N, {0, 0, 0});

    Solver solver {M, m, T/N, N_small, N, pos_init, vel_init, big_positions, big_velocities};

    auto t0 = std::chrono::steady_clock::now();
    solver.solve();
    auto t1 = std::chrono::steady_clock::now();
    std::cout << "Time: " << std::chrono::duration<double>{t1-t0}.count() << "s" << std::endl;

    std::ofstream file;
    file.open(filename);
    file.precision(4);
    file << std::left;
    for (unsigned i = 0; i < N; i++) {
        file.width(5);
        file << i*T/N;
        for (int j : {0, 1, 2}) {
            file.width(14);
            file << big_positions[i][j];
        }
        for (int j : {0, 1, 2}) {
            if (j < 2) {
                file.width(14);
            }
            file << big_velocities[i][j];
        }
        file << "\n";
    }
    file.close();

    return 0;
}