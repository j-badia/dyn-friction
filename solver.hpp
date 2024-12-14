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
    void calc_accelerations(const vec& positions);
    void multiply(vec& v, double a);
    void add_to_first(vec& v1, const vec& v2);
    void add_to_first_then_mult(vec& v1, const vec& v2, const double a);
    void mult_then_add_to_first(vec& v1, const vec& v2, const double a);

    public:

    Solver(double _M, double _m, double _step, unsigned _N_small, unsigned _N_steps, const vec& init_pos, const vec& init_vel, vec3& _big_pos, vec3& _big_vel);
    void solve();
};