
#include "utils.cpp"
#include "math.h"
#include <vector>
#include <iostream>

double F(double m) {
    return m;
}

double F_deriv(double m) {
    return 1.0;
}

double G(double m) {
    return m * (1.0 - m);
}

double G_deriv(double m) {
    return 1.0 - 2.0 * m;
}

double E(double m) {
    return m;
}

double E_deriv(double m) {
    return 1.0;
}

double initial_density(double x) {
    return 0.95 * exp(-20.0 * x * x);
}

// double initial_density(double x) {
//     return 0.03;
// }
//
/*
double I_T(std::vector<Vector> ms, std::vector<Vector> vs, double dt, double dx) {
    double result = 0.0;
    for (int time = 0; time < ms.size(); time++) {
        for (int x = 0; x < ms[0].len; x++) {
            result += (F_deriv(ms[time][x]) * abs(vs[time][x]) * abs(vs[time][x])
                + E(ms[time][x])) * dx * dt * 0.5;
        }
    }
    return result;
}
*/

double I_T_diff(std::vector<Vector> ms, std::vector<Vector> vs, double dt, double dx) {
	double IT1 = 0.0; double IT2 = 0.0;
		for(int x = 0; x < ms[0].len; x++) {
			IT1 += (F_deriv(ms[time + 1][x] * abs(vs[time + 1][x]) *
			abs(vs[time + 1][x] + E(ms[time + 1][x])) *dx * dt * 0.5);
			IT2 += (F_deriv(ms[time][x] * abs(vs[time][x]) *
			abs[time][x] + E(ms[time][x])) * dx * dt * 0.5);
		}
	return abs(IT1 - IT2);
}

int main(int argc, char **argv) {
    double d = 0.005;
    double p = 1.0;
    double T = 2.0;

    int mspace = 20;
    int ntime = 20;
    int niters = 20;

    double dx = 2.0 / ((double) mspace);
    double dt = T / ((double) ntime);

    double k = 0.01;

    double e1 = 0.0001;

    std::vector<Vector> ms;
    for (int i = 0; i < ntime; i++) {
        ms.push_back(Vector(mspace));
    }

    // Set initial density
    for (int x = 0; x < mspace; x++) {
        ms[0][x] = initial_density(((double) x) / (double) mspace);
    }

    std::vector<Vector> us;
    for (int i = 0; i < ntime; i++) {
        us.push_back(Vector(mspace));
    }

    std::vector<Vector> vs;
    for (int i = 0; i < ntime; i++) {
        vs.push_back(Vector(mspace));
    }
    
    for (int time = 0; time < ntime; time++) {
        for (int x = 0; x < mspace; x++) {
            ms[time][x] = initial_density(((double) x) / (double) mspace);
        }
    }

    double I_T_last;

    bool first = true;

    for (int j = 0; j < niters; j++) {
        for (int time = 0; time < ntime - 2; time++) {
            int r = mspace - 1;
            // Assign the velocities depending on the adjoint variables
            double nabla_u;
            Vector vs_old (vs[time]);
            for (int x = 1; x < mspace - 1; x++) {
                nabla_u = (us[time + 1][x + 1] - us[time + 1][x]) / (dx);
                vs[time + 1][x] = G(ms[time + 1][x]) / F(ms[time + 1][x]) * nabla_u;
            }
            nabla_u = (us[time + 1][1] - us[time + 1][0]) / dx;
            vs[time + 1][0] = G(ms[time + 1][0]) / F(ms[time + 1][0]) * nabla_u;
            nabla_u = (us[time + 1][r] - us[time + 1][r - 1]) / dx;
            vs[time + 1][r] = G(ms[time + 1][r]) / F(ms[time + 1][r]) * nabla_u;
 
/*~~~~~~~~~~~~~~~~~~~~~*/  	    
	    // Step 1
            // eq 19a
            for (int x = 1; x < mspace - 1; x++) {
                double delta_m =
                    (ms[time + 1][x + 2] - 2* ms[time + 1][x + 1] + ms[time + 1][x]) / dx / dx;
                //double other_term = ((F(ms[time][x + 1]) * vs[time][x + 1]) - ((F(ms[time][x - 1]) * vs[time][x - 1]))) / (2.0 * dx);
                // double delt_m = other_term - delta_m;
		double other_term = 0.0; //because of the galerkin weights
                ms[time][x] = ms[time + 1][x] - (d * delta_m - other_term) * dt;
            }
            // eq 19b
            // Newton's method
            // Left boundary
	    // //TODO: double check this
            double m0 = ms[time][0];
            const double m1 = ms[time][1];
            const double v0 = vs[time][0];
            const double v1 = vs[time][1];
            for (int i = 0; i < 3; i++) {
                double N = -d * (m1 - m0) / dx + v0 *(F(m1) - F(m0)) / dx - p * m0 * v0;
                double N_prime = d / dx - v0 * F_deriv(m0) / dx - p * v0;
                m0 -= N / N_prime;
                // printf("\n\nm0: should converge: %f, %f\n\n\n", m0, N_prime);
            }
            ms[time + 1][0] = m0;
            // Right boundary
            double mr = ms[time][r];
            const double mrm1 = ms[time][r - 1];
            const double vr = vs[time][r];
            const double vrm1 = vs[time][r - 1];
            for (int i = 0; i < 3; i++) { 
                double N = -d * (mr - mrm1) / dt + vr * (F(mr) - F(mrm1)) / dx - p * mrm1 * vr;
                double N_prime = d / dx - vr * F_deriv(mrm1) / dx - p * vr;
                mr -= N / N_prime;
            }
            ms[time + 1][r] = mr;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~*/

            // Step 2
            for (int x = 1; x < mspace - 1; x++) {
                double delta_u =
                    (us[time + 1][x] - 2 * us[time + 1][x - 1] + us[time + 1][x -2]) / dx / dx;
                double nabla_u = (us[time + 1][x] - us[time + 1][x - 1]) / (dx);
                double m = ms[time + 1][x];
                double v = vs[time + 1][x];
                us[time][x] = us[time + 1][x] - 
                    (0.5 * (F_deriv(m) * abs(v) * abs(v) + E_deriv(m)) -
                    (d * delta_u + G_deriv(m) * v * nabla_u)) * dt;
            }
            // Left boundary
            // no newton's method!
	    // TODO: come back to  
            us[time + 1][0] = ((-d * us[time + 1][1]) / dx - p * vs[time + 1][0] / dx);
            us[time][r] = ((d * us[time][r - 1]) / dx) / (p * vs[time][r] + d / dx);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/	    
            // Step 3
            for (int x = 1; x < mspace - 1; x++) {
                nabla_u = (us[time + 2][x + 1] - us[time + 1][x]) * (dx);
                vs[time + 1][x] = vs_old[x] - k * (F(ms[time + 1][x]) * vs[time + 1][x] - G(ms[time + 1][x]) * nabla_u);
            }
	    
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

            // Step 4, check for convergence
	    // Want the [time+1] values to be I_T_last and the [time] values to be I_T_new
/*
            if (first) {
                first = false;
                I_T_last = I_T(ms, vs, dt, dx);
	   
            }else{ 
           	double I_T_new = I_T(ms, vs, dt, dx);
            	double error = abs(I_T_last - I_T_new);
            	I_T_last = I_T_new;
*/
	        double error = IT_diff(ms, vs, dt, dx);
	
		fprintf(stderr, "Error: %2.15e\n", error);
            	if (error < e1) {
                	fprintf(stderr, "Converged!\n");
                	goto END;
            	}
            	if (error > 1e30) {
                	fprintf(stderr, "Diverged!\n");
                	goto END;
           	 }
	    }
        }
    }

END:

    // Print results
    printf("ms:\n");
    for (int time = 0; time < ntime; time++) {
        printf("%d: ", time);
        for (int x = 0; x < mspace; x++) {
            printf("%f, ", ms[time][x]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("vs:\n");
    for (int time = 0; time < ntime; time++) {
        printf("%d: ", time);
        for (int x = 0; x < mspace; x++) {
            printf("%f, ", vs[time][x]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("us:\n");
    for (int time = 0; time < ntime; time++) {
        printf("%d: ", time);
        for (int x = 0; x < mspace; x++) {
            printf("%f, ", us[time][x]);
        }
        printf("\n");
    }
    printf("\n");
}
