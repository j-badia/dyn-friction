#include <cassert>
#include <cmath>
#include <array>
#include <vector>
#include <iostream>

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

    constexpr static double G = 4.3e-3; // In pc (km/s)^2 / M_sun
    double M; // Mass of big particle
    double m; // Mass of small particles
    double step; // Time step size
    unsigned N_small; // Number of small particles
    unsigned curr_step {0};
    
    //These vectors hold the position and velocities of all the particles. The big mass comes first, so the size is 3*(N_small+1).
    vec curr_positions;
    vec curr_velocities;
    vec next_positions;
    vec next_velocities;
    vec accels;

    // These Nx3 vectors hold the data for the big particle. The caller must preallocate but doesn't need to put the initial values in.
    vec3& big_positions;
    vec3& big_velocities;
    
    // Will calculate the accelerations by taking into account only the forces between the first particle (with mass M) and all the others (with mass m), and store them in accels.
    void calc_accelerations(const vec& positions);
    void multiply(vec& v, double a);
    void add_to_first(vec& v1, const vec& v2);
    void add_to_first_then_mult(vec& v1, const vec& v2, const double a);
    void mult_then_add_to_first(vec& v1, const vec& v2, const double a);

    public:

    // init_pos is 3(N_small+1), big mass first
    Solver(double _M, double _m, double _step, const vec& init_pos, const vec& init_vel, vec3& _big_pos, vec3& _big_vel);
    void solve(unsigned N_steps);
    void next_step();
};