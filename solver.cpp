#include "solver.hpp"

void Solver::calc_accelerations(const vec& positions) {
    accels.assign(accels.size(), 0);
    for (unsigned i = 1; i < N_small+1; i++) {
        std::array<double, 3> pos_rel {
            positions[3*i] - positions[0],
            positions[3*i+1] - positions[1],
            positions[3*i+2] - positions[2]
        };
        double d2 = pos_rel[0]*pos_rel[0]+pos_rel[1]*pos_rel[1]+pos_rel[2]*pos_rel[2];
        double inv_d3 = G / std::sqrt(d2*d2*d2);
        for (int j : {0, 1, 2}) {
            accels[j] += m * inv_d3 * pos_rel[j];
            accels[3*i+j] -= M * inv_d3 * pos_rel[j];
        }
    }
}

void Solver::multiply(vec& v, double a) {
    for (auto& it : v) {
        it *= a;
    }
}

void Solver::add_to_first(vec& v1, const vec& v2) {
    unsigned N = v1.size();
    for (unsigned i = 0; i < N; i++) {
        v1[i] += v2[i];
    }
}

void Solver::add_to_first_then_mult(vec& v1, const vec& v2, const double a) {
    unsigned N = v1.size();
    for (unsigned i = 0; i < N; i++) {
        v1[i] += v2[i];
        v1[i] *= a;
    }
}

void Solver::mult_then_add_to_first(vec& v1, const vec& v2, const double a) {
    unsigned N = v1.size();
    for (unsigned i = 0; i < N; i++) {
        v1[i] += a * v2[i];
    }
}

Solver::Solver(double _M, double _m, double _step, const vec& init_pos, const vec& init_vel, vec3& _big_pos, vec3& _big_vel)
    : M{_M}, m{_m}, step{_step}, big_positions{_big_pos}, big_velocities{_big_vel}
{
    assert(init_pos.size() == init_vel.size());
    //assert(big_positions.size() == big_velocities.size()); Is this correct if they're reserved but not initialized?
    N_small = init_pos.size()/3 - 1;
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

void Solver::next_step() {
    next_positions = curr_positions;
    next_velocities = curr_velocities;

    mult_then_add_to_first(next_positions, curr_velocities, step/2);
    calc_accelerations(next_positions);
    mult_then_add_to_first(next_velocities, accels, step);
    //next_positions = curr_positions;
    mult_then_add_to_first(next_positions, next_velocities, step/2);

    curr_step++;
    for (int j : {0, 1, 2}) {
        big_positions[curr_step][j] = next_positions[j];
        big_velocities[curr_step][j] = next_velocities[j];
    }

    std::swap(curr_positions, next_positions);
    std::swap(curr_velocities, next_velocities);
}

void Solver::solve(unsigned N_steps) {
    unsigned print_interval = (unsigned) N_steps/100;
    for (unsigned i = 1; i < N_steps; i++) {
        if (i % print_interval == 0) {
            std::cout << "\rt = " << i*step << std::flush;
        }
        next_step();
    }
    std::cout << "\n";
}