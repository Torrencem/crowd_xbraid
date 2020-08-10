#include <math.h>

double nice_curve(double x, double c) {
    double exponent = 1.0 / (1.0 + (x - c) * (x - c));
    double pow_of_2 = powf(2.0, exponent - 1.0);
    return powf(pow_of_2, 80.0);
}

/// Initial condition for rho
double initial_condition(double x) {
    return 0.8 * nice_curve(x, 0.3) + 0.3 * nice_curve(x, 0.7);
}

/// Final condition for rho
double final_condition(double x) {
    return 0.1 * nice_curve(x, 0.3) + 0.8 * nice_curve(x, 0.7);
}
