#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include <cstdint>

const uint64_t N {10000};

using namespace std::chrono;

/* 
 *  Testing ways of calculating x^(-3/2)
 */


//The double loop is because for some reason the loop was being skipped if N was large
void pow_sqrt(double x) {
    double total = 0;
    for (uint64_t j = 0; j < N; j++) {
        if (j == 0) {
            std::cout << "pow_sqrt in loop" << "\n";
        }
        for (uint64_t i = 0; i < N; i++) {
            total += 1.0 / pow(sqrt(x), 3);
        }
    }
}

void manual_cube_sqrt(double x) {
    double total = 0;
    for (uint64_t j = 0; j < N; j++) {
        if (j == 0) {
            std::cout << "manual_cube_sqrt in loop" << "\n";
        }
        for (uint64_t i = 0; i < N; i++) {
            double y = x*x*x;
            total += 1.0 / sqrt(y);
        }
    }
}

void only_pow(double x) {
    double total = 0;
    for (uint64_t j = 0; j < N; j++) {
        if (j == 0) {
            std::cout << "only_pow in loop" << "\n";
        }
        for (uint64_t i = 0; i < N; i++) {
            total += pow(x, -1.5);
        }
    }
}

int main(int argc, char* argv[]) {
    double x = 20;

    auto t0 = steady_clock::now();
    pow_sqrt(x);
    auto t1 = steady_clock::now();
    std::cout << "Pow and sqrt functions: " << duration<double>{t1-t0}.count() << "s\n";

    t0 = steady_clock::now();
    manual_cube_sqrt(x);
    t1 = steady_clock::now();
    std::cout << "Manual cube plus sqrt: " << duration<double>{t1-t0}.count() << "s\n";

    t0 = steady_clock::now();
    only_pow(x);
    t1 = steady_clock::now();
    std::cout << "Pow with -1.5 exponent: " << duration<double>{t1-t0}.count() << "s\n";

    return 0;
}