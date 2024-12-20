/*
    Units:
    distance in parsecs
    velocity in km/s
    time in pc/(km/s) = 978,462 yr
    mass in solar masses

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
#include <sstream>
#include <map>
#include <chrono>

#include "solver.hpp"

const std::string default_params_file {"params"};

using smap = std::map<std::string, std::string>;

smap read_config(std::string fname) {
    smap config;
    std::ifstream file;
    file.open(fname);
    std::string name, value;
    while (file.peek() != EOF) {
        if (file.peek() == '#') {
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        std::getline(file, name, '=');
        std::getline(file, value);
        config[name] = value;
    }
    return config;
}

/*
 *  Reads a file consisting of two lines of single space separated numbers, and
 *  converts them into two vectors.
 */
std::pair<vec, vec> read_initial(std::string fname) {
    std::ifstream file;
    file.open(fname);
    std::string pos_s, vel_s;
    std::getline(file, pos_s);
    std::getline(file, vel_s);
    file.close();

    vec pos, vel;
    std::stringstream ss {pos_s};
    std::string elem;
    while (std::getline(ss, elem, ' ')) {
        pos.emplace_back(std::stod(elem));
    }
    ss = std::stringstream(vel_s);
    while (std::getline(ss, elem, ' ')) {
        vel.emplace_back(std::stod(elem));
    }
    
    return std::make_pair(pos, vel);
}

/*
args file:
initial: initial conditions file, two lines: one for position and one for velocity
output: output file to store trajectory of big mass
M: big mass
m: small mass
T: integration time interval
dt: integration time step
*/
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: df.exe args-file\n";
        return 1;
    }
    smap config = read_config(argv[1]);
    for (const std::string& key : {"initial", "output", "M", "m", "T", "dt"}) {
        if (config.count(key) == 0) {
            std::cout << "Missing element " << key << " in args.\n";
            return 1;
        }
    }

    double M, m, T, dt;
    M = std::stod(config["M"]);
    m = std::stod(config["m"]);
    T = std::stod(config["T"]);
    dt = std::stod(config["dt"]);
    
    auto [init_pos, init_vel] = read_initial(config["initial"]);
    unsigned N_small = init_pos.size();
    if (N_small != init_vel.size()) {
        std::cout << "Initial position and velocity have different lengths.\n";
        return 1;
    }

    unsigned N_steps = (unsigned) (T/dt);
    vec3 big_pos(N_steps+1, {0, 0, 0});
    vec3 big_vel(N_steps+1, {0, 0, 0});
    Solver solver {M, m, dt, init_pos, init_vel, big_pos, big_vel};

    auto t0 = std::chrono::steady_clock::now();

    unsigned print_interval = N_steps/100;
    for (unsigned i = 1; i < N_steps+1; i++) {
        if (i % print_interval == 0) {
            std::cout << "\rt = " << i*dt << std::flush;
        }
        solver.next_step();
    }
    std::cout << "\n";

    auto t1 = std::chrono::steady_clock::now();
    std::cout << "Time: " << std::chrono::duration<double>{t1-t0}.count() << "s" << std::endl;

    std::ofstream file;
    file.open(config["output"]);
    file.precision(4);
    file << std::left;
    for (unsigned i = 0; i < N_steps+1; i++) {
        file.width(5);
        file << i*dt;
        for (int j : {0, 1, 2}) {
            file.width(14);
            file << big_pos[i][j];
        }
        for (int j : {0, 1, 2}) {
            if (j < 2) {
                file.width(14);
            }
            file << big_vel[i][j];
        }
        file << "\n";
    }
    file.close();

    return 0;
}

int old_main(int argc, char* argv[]) {
    double M, m, distance, v0, d0, T;
    unsigned N_width, N_length, N;
    std::string filename;
    smap config;
    if (argc == 1 || argc == 2) {
        std::string params_file;
        if (argc == 1) {
            params_file = default_params_file;
        } else if (argc == 2) {
            params_file = argv[1];
        }
        config = read_config(params_file);
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

    Solver solver {M, m, T/N, pos_init, vel_init, big_positions, big_velocities};

    auto t0 = std::chrono::steady_clock::now();
    solver.solve(N);
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