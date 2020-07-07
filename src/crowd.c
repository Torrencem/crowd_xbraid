/*BHEADER**********************************************************************
 * Written by Isaiah Meyers, Joseph Munar, Eric Neville, Tom Overman
 *
 * This file is part of XBraid. For support, post issues to the XBraid Github
 *page.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free
 *Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 *ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 *FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 *Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 *along with this program; if not, write to the Free Software Foundation, Inc.,
 *59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/

/**
 * Example:       model_problem.c
 *
 * Interface:     C
 *
 * Requires:      Lapacke
 *
 * Compile with:  make model_problem
 *
 * Description:   Solves a crowd control problem in time-parallel:
 *
 *                min   1/2 \int_0^T \int_0^1 |w(x, t)|^2+v(x, t)^2 dx dt
 *
 *                s.t.  d/dt rho = sigma^2/2 d^2/dx^2 rho - d/dx K(rho)*w
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "braid.h"
#include "braid_test.h"
#define PI 3.14159265

#include "lapacke.h"
#include "utils.c"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct {
    int myid;   /* Rank of the processor */
    int ntime;  /* Total number of time-steps (starting at time 0) */
    int mspace; /* Number of space points included in our state vector */
                /* So including boundaries we have M+2 space points */
    double dx;  /* Spatial mesh spacing */

    double sigma2; /* Sigma^2 parameter for crowd simulation */
    double k_coeff; /* We assume K(rho) = k_coeff * rho + k_const */
    double k_const;

    int ilower;  /* Lower index for my proc */
    int iupper;  /* Upper index for my proc */
    int npoints; /* Number of time points on my proc */
    FILE *rhofile, *wfile, *zfile;
} my_App;

/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct {
    Vector rho;
    Vector w;
    Vector z; /* Lagrange multiplier weights */
} my_Vector;

double initial_condition_rho(const double x) {
	// \frac{\left(2^{\frac{1}{1\ +\ \left(4x-\frac{4T}{2}\right)^{2}}}\ -\ 1\right)}{1.5}
	// (2^(1 / (1 + (4x - (4T)/2))) - 1) / 1.5
	// https://www.desmos.com/calculator/crvh8cl2ta
	// Endpoint in space. TODO this should be probably #defined
	double T = 1.0;
	return (pow(2.0, 1.0 / (1.0 + (4.0 * x - (4.0 * T) / 2.0))) - 1.0) / 1.5;
}

double initial_condition_w(const double x) {
	// sqrt(F(p))v, but since v is initially 0, w is initially 0
	return 0.0;
}

double get_rho(const my_App *app, const Vector rho, const int n){
	if(0 <= n && n < app->mspace){
		return rho[n];
	} else {
		return 0.0; //TODO: Boundary conditions
	}
}

double get_w(const my_App *app, const Vector w, const int n){
	if(0 <= n && n < app->mspace){
		return w[n];
	} else {
		return 0.0; //TODO: Boundary conditions
	}
}

#define K(r) (app->k_coeff * r + app->k_const)
Vector apply_Phi(const my_App *app, const Vector rho, const Vector w,
                 const double t, const double dt) {

    double gamma = dt/(2.0 * app->dx);
    Vector rho_new = zero_vector(app->mspace);
    
    for (int j = 0; j < app->mspace; j++){
    	rho_new[j] = (app->sigma2 / app->dx) * (get_rho(app, rho, j+1) - 2.0*get_rho(app, rho, j) + get_rho(app, rho, j-1));
        rho_new[j] -= K(get_rho(app, rho, j)) * (get_w(app, w, j+1) - get_w(app, w, j-1));
        rho_new[j] += get_w(app, w, j) * (K(get_rho(app, rho, j+1)) - K(get_rho(app, rho, j-1)));
        rho_new[j] *= gamma;
        rho_new[j] += get_rho(app, rho, j);
    }

    return rho_new;
}


void compute_R_matrix(const my_App *app, const Vector rho, const int n,
                      Vector *LL, Vector *L, Vector *LU, const double dt,
                      const double t) {

    double coeff = -dt/(2.0 * app->dx);

    *LL = zero_vector(n - 1);
    for (int i = 0; i < n - 1; i++) {
        (*LL)[i] = coeff * K(get_rho(app, rho, i+1));
    }

    *L = zero_vector(n);
    for (int i = 0; i < n; i++) {
        (*L)[i] = coeff * (K(get_rho(app, rho, i+1)) - K(get_rho(app, rho, i-1)));
    }

    *LU = zero_vector(n - 1);
    for (int i = 0; i < n - 1; i++) {
        (*LU)[i] = coeff * K(get_rho(app, rho, i));
    }
}


void compute_S_matrix(const my_App *app, const Vector rho, const Vector w, const int n,
                      Vector *LL, Vector *L, Vector *LU, const double dt,
                      const double t) {

    double s2dx = app->sigma2 / app->dx;
    double coeff = dt/(2.0 * app->dx);

    *LL = zero_vector(n - 1);
    for (int i = 0; i < n - 1; i++) {
        (*LL)[i] = coeff * (s2dx + K(get_rho(app, rho, i)) * w[i+1]);
    }

    *L = zero_vector(n);
    for (int i = 0; i < n; i++) {
        (*L)[i] = 1.0 + coeff * (-2.0 * s2dx + K(get_rho(app, rho, i)) * (get_w(app, w, i+1) - get_w(app, w, i-1)));
    }

    *LU = zero_vector(n - 1);
    for (int i = 0; i < n - 1; i++) {
        (*LU)[i] = coeff * (s2dx - K(get_rho(app, rho, i+1)) * w[i]);
    }
}

/*------------------------------------*/

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */

int my_TriResidual(braid_App app, braid_Vector uleft, braid_Vector uright,
                   braid_Vector f, braid_Vector r, braid_TriStatus status) {
    double t, tprev, tnext, dt;
    int level, index;

    braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
    braid_TriStatusGetLevel(status, &level);
    braid_TriStatusGetTIndex(status, &index);

    /* Get the time-step size */
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }

    Vector rho_res;
    // Compute u
    if (uleft == NULL) {
        // Collect initial conditions
        Vector rho = zero_vector(app->mspace);
        Vector w = zero_vector(app->mspace);
        for (int x = 0; x < app->mspace; x++) {
            rho[x] = initial_condition_rho(x);
            w[x] = initial_condition_w(x);
        }
        rho_res = apply_Phi(app, rho, w, t, dt);
    } else {
        rho_res = apply_Phi(app, uleft->rho, uleft->w, t, dt);
    }
    vec_axpy(app->mspace, 1.0, r->rho, -1.0, rho_res);

    Vector w_res = zero_vector(app->mspace);
    if (uright != NULL) {
        vec_copy(app->mspace, uright->w, v_res);
        Vector LL = NULL;
        Vector L = NULL;
        Vector LU = NULL;
        compute_R_matrix(app, r->rho, app->mspace, &LL, &L, &LU, dt, t);
        matmul_tridiag(LL, L, LU, app->mspace, &v_res);
    }
    vec_axpy(app->mspace, app->dx * dt, r->w, -1.0, v_res);

    Vector z_res = zero_vector(app->mspace);
    if (uright != NULL) {
        vec_copy(app->mspace, uright->w, w_res);
        Vector LL = NULL;
        Vector L = NULL;
        Vector LU = NULL;
        compute_S_matrix(app, r->rho, r->w, app->mspace, &LL, &L, &LU, dt, t);
        matmul_tridiag(LL, L, LU, app->mspace, &w_res);
    }
    vec_axpy(app->mspace, -1.0 * app->dx * dt, r->rho, 1.0, w_res);

    if (f != NULL){
        vec_axpy(app->mspace, -1.0, f->rho, 1.0, rho_res);
        vec_axpy(app->mspace, -1.0, f->w, 1.0, w_res);
        vec_axpy(app->mspace, -1.0, f->z, 1.0, z_res);
    }

    vec_copy(app->mspace, rho_res, r->rho);
    vec_copy(app->mspace, w_res, r->w);
    vec_copy(app->mspace, z_res, r->z);

    return 0;
}

/*------------------------------------*/

/* Solve A(u) = f */

int my_TriSolve(braid_App app, braid_Vector uleft, braid_Vector uright,
                braid_Vector fleft, braid_Vector fright, braid_Vector f,
                braid_Vector u, braid_TriStatus status) {
    double dx = (app->dx);

    double t, tprev, tnext, dt;

    /* Get the time-step size */
    braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }

    Vector rho_new;
    if (uleft == NULL) {
        // Collect initial conditions
        Vector rho = zero_vector(app->mspace);
        Vector w = zero_vector(app->mspace);
        for (int x = 0; x < app->mspace; x++) {
            rho[x] = initial_condition_w(x);
            w[x] = initial_condition_rho(x);
        }
        u_new = apply_Phi(app, rho, w, t, dt);
    } else {
        rho_new = apply_Phi(app, uleft->rho, uleft->w, t, dt);
    }

    Vector w_new = zero_vector(app->mspace);
    if (uright != NULL) {
        Vector LL = NULL;
        Vector L = NULL;
        Vector LU = NULL;
        compute_R_matrix(app, rho_new, app->mspace, &LL, &L, &LU, dt, t);
        vec_copy(app->mspace, uright->z, w_new);
        matmul_tridiag(LL, L, LU, app->mspace, &w_new);
        vec_scale(app->mspace, 1.0 / (dx * dt), w_new);
    }
    
    Vector z_new = zero_vector(app->mspace);
    if (uright != NULL) {
        vec_copy(app->mspace, uright->w, w_new);
        Vector LL = NULL;
        Vector L = NULL;
        Vector LU = NULL;
        compute_S_matrix(app, rho_new, w_new, app->mspace, &LL, &L, &LU, dt, t);
        matmul_tridiag(LL, L, LU, app->mspace, &w_new);
    }
    vec_axpy(app->mspace, -1.0 * dx * dt, u_new, 1.0, w_new);

    vec_copy(app->mspace, rho_new, u->rho);
    vec_copy(app->mspace, w_new, u->w);
    vec_copy(app->mspace, z_new, u->z);

    /* no refinement */
    braid_TriStatusSetRFactor(status, 1);

    return 0;
}

/*------------------------------------*/

/* This is only called from level 0 */

int my_Init(braid_App app, double t, braid_Vector *u_ptr) {
    int i;
    my_Vector *u;
    int mspace = (app->mspace);

    /* Allocate the vector */
    u = (my_Vector *)malloc(sizeof(my_Vector));
    u->rho = zero_vector(mspace);
    u->w = zero_vector(mspace);
    u->z = zero_vector(mspace);

    for (i = 0; i <= mspace - 1; i++) {
        u->rho[i] = ((double)braid_Rand()) / braid_RAND_MAX;
        u->w[i] = ((double)braid_Rand()) / braid_RAND_MAX;
        u->z[i] = ((double)braid_Rand()) / braid_RAND_MAX;
    }

    *u_ptr = u;

    return 0;
}

/*------------------------------------*/

int my_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr) {
    int mspace = (app->mspace);
    my_Vector *v;

    /* Allocate the vector */
    v = (my_Vector *)malloc(sizeof(my_Vector));

    v->rho = zero_vector(app->mspace);
    v->w = zero_vector(app->mspace);
    v->z = zero_vector(app->mspace);

    vec_copy(mspace, (u->rho), (v->rho));
    vec_copy(mspace, (u->w), (v->w));
    vec_copy(mspace, (u->z), (v->z));
    *v_ptr = v;

    return 0;
}

/*------------------------------------*/

int my_Free(braid_App app, braid_Vector u) {
    free(u->rho);
    free(u->w);
    free(u->z);
    free(u);

    return 0;
}

/*------------------------------------*/

int my_Sum(braid_App app, double alpha, braid_Vector x, double beta,
           braid_Vector y) {

    vec_axpy((app->mspace), alpha, (x->rho), beta, (y->rho));
    vec_axpy((app->mspace), alpha, (x->w), beta, (y->w));
    vec_axpy((app->mspace), alpha, (x->z), beta, (y->z));
    return 0;
}

/*------------------------------------*/

int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr) {
    int i;
    double dot = 0.0;
    int mspace = (app->mspace);
    for (i = 0; i < mspace; i++) {
        dot += (u->rho)[i] * (u->rho)[i];
        dot += (u->w)[i] * (u->w)[i];
        dot += (u->z)[i] * (u->z)[i];
    }
    *norm_ptr = sqrt(dot);

    return 0;
}

/*------------------------------------*/

int my_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus) {
    int done, index;
    int mspace = (app->mspace);

    /* Print solution to file if simulation is over */
    braid_AccessStatusGetDone(astatus, &done);
    if (done) {
        int j;

        braid_AccessStatusGetTIndex(astatus, &index);

        fprintf((app->rhofile), "%05d: ", index);
        fprintf((app->wfile), "%05d: ", index);
        fprintf((app->zfile), "%05d: ", index);
        for (j = 0; j < (mspace - 1); j++) {
            fprintf((app->rhofile), "% 1.14e, ", (u->rho[j]));
            fprintf((app->wfile), "% 1.14e, ", (u->w[j]));
            fprintf((app->zfile), "% 1.14e, ", (u->z[j]));
        }
        fprintf((app->rhofile), "% 1.14e\n", (u->rho[j]));
        fprintf((app->wfile), "% 1.14e\n", (u->w[j]));
        fprintf((app->zfile), "% 1.14e\n", (u->z[j]));
    }

    return 0;
}

/*------------------------------------*/

int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus bstatus) {
    *size_ptr = 3 * (app->mspace) * sizeof(double);
    return 0;
}

/*------------------------------------*/

int my_BufPack(braid_App app, braid_Vector u, void *buffer,
               braid_BufferStatus bstatus) {
    double *dbuffer = buffer;

    vec_copy((app->mspace), (u->rho), dbuffer);
    dbuffer += (app->mspace);
    vec_copy((app->mspace), (u->w), dbuffer);
    dbuffer += (app->mspace);
    vec_copy((app->mspace), (u->z), dbuffer);
    braid_BufferStatusSetSize(bstatus, 3 * (app->mspace) * sizeof(double));

    return 0;
}

/*------------------------------------*/

int my_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr,
                 braid_BufferStatus bstatus) {
    my_Vector *u = NULL;
    double *dbuffer = buffer;

    /* Allocate memory */
    u = (my_Vector *)malloc(sizeof(my_Vector));
    vec_create((app->mspace), &(u->rho));
    vec_create((app->mspace), &(u->w));
    vec_create((app->mspace), &(u->z));

    /* Unpack the buffer */
    vec_copy((app->mspace), dbuffer, (u->rho));
    dbuffer += (app->mspace);
    vec_copy((app->mspace), dbuffer, (u->w));
    dbuffer += (app->mspace);
    vec_copy((app->mspace), dbuffer, (u->z));

    *u_ptr = u;
    return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    braid_Core core;
    my_App *app;

    double tstart, tstop, dt, dx, start, end;
    int rank, ntime, mspace, arg_index;
    int max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
    int access_level, print_level;
    double tol;
    double time;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Define space domain. Space domain is between 0 and 1, mspace defines the
     * number of steps */
    mspace = 16;
    ntime = 16;

    /* Define some Braid parameters */
    max_levels = 2;
    min_coarse = 1;
    nrelax = 1;
    nrelaxc = 10;
    maxiter = 50;
    cfactor = 2;
    tol = 1.0e-6;
    access_level = 2;
    print_level = 2;

    /* Parse command line */
    arg_index = 1;
    while (arg_index < argc) {
        if (strcmp(argv[arg_index], "-help") == 0) {
            printf("\n");
            printf(" Solves the advection-diffusion model problem \n\n");
            printf("  min  1/2 \\int_0^T\\int_0^1 (u(x,t)-ubar(x))^2 + "
                   "alpha*v(x,t)^2  dxdt \n\n");
            printf("  s.t.  u_t + u_x - nu*u_xx = v(x,t) \n");
            printf("        u(0,t) = u(1,t) = 0 \n\n");
            printf("        u(x,0) = u0(x) \n");
            printf("  -tstop <tstop>          : Upper integration limit for "
                   "time\n");
            printf("  -ntime <ntime>          : Num points in time\n");
            printf("  -mspace <mspace>        : Num points in space\n");
            printf("  -ml <max_levels>        : Max number of braid levels \n");
            printf("  -num  <nrelax>          : Num F-C relaxations\n");
            printf("  -nuc <nrelaxc>          : Num F-C relaxations on "
                   "coarsest grid\n");
            printf("  -mi <maxiter>           : Max iterations \n");
            printf("  -cf <cfactor>           : Coarsening factor \n");
            printf("  -tol <tol>              : Stopping tolerance \n");
            printf("  -access <access_level>  : Braid access level \n");
            printf("  -print <print_level>    : Braid print level \n");
            printf("\n");
            exit(1);
        } else if (strcmp(argv[arg_index], "-ntime") == 0) {
            arg_index++;
            ntime = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-tstop") == 0) {
            arg_index++;
            tstop = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-mspace") == 0) {
            arg_index++;
            mspace = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-ml") == 0) {
            arg_index++;
            max_levels = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-num") == 0) {
            arg_index++;
            nrelax = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-nuc") == 0) {
            arg_index++;
            nrelaxc = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-mi") == 0) {
            arg_index++;
            maxiter = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-cf") == 0) {
            arg_index++;
            cfactor = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-tol") == 0) {
            arg_index++;
            tol = atof(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-access") == 0) {
            arg_index++;
            access_level = atoi(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-print") == 0) {
            arg_index++;
            print_level = atoi(argv[arg_index++]);
        } else {
            printf("ABORTING: incorrect command line parameter %s\n",
                   argv[arg_index]);
            return (0);
        }
    }

    /* Define the space step */
    dx = (double)1 / (mspace + 1);

    /* Define time domain and step */
    tstart = 0.0; /* Beginning of time domain */
    tstop = 1.0;  /* End of time domain*/
    dt = (tstop - tstart) / ntime;

    /* Set up the app structure */
    app = (my_App *)malloc(sizeof(my_App));
    app->myid = rank;
    app->ntime = ntime;
    app->mspace = mspace;
    app->dx = dx;

    char filename[255];
    sprintf(filename, "%s.%03d", "modelproblem.out.rho", (app->myid));
    (app->rhofile) = fopen(filename, "w");
    sprintf(filename, "%s.%03d", "modelproblem.out.w", (app->myid));
    (app->wfile) = fopen(filename, "w");
    sprintf(filename, "%s.%03d", "modelproblem.out.z", (app->myid));
    (app->zfile) = fopen(filename, "w");

    /* Initialize XBraid */
    braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, dt, tstop, ntime - 1,
                       app, my_TriResidual, my_TriSolve, my_Init, my_Clone,
                       my_Free, my_Sum, my_SpatialNorm, my_Access, my_BufSize,
                       my_BufPack, my_BufUnpack, &core);

    /* Set some XBraid(_Adjoint) parameters */
    braid_SetMaxLevels(core, max_levels);
    braid_SetMinCoarse(core, min_coarse);
    braid_SetNRelax(core, -1, nrelax);
    if (max_levels > 1) {
        braid_SetNRelax(core, max_levels - 1,
                        nrelaxc); /* nrelax on coarsest level */
    }
    braid_SetCFactor(core, -1, cfactor);
    braid_SetAccessLevel(core, access_level);
    braid_SetPrintLevel(core, print_level);
    braid_SetMaxIter(core, maxiter);
    braid_SetAbsTol(core, tol);

    /* Parallel-in-time TriMGRIT simulation */
    start = clock();
    braid_Drive(core);
    end = clock();

    /* Print runtime to file (for runtime comparisons)*/
    time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total Run Time: %f s \n", time);

    free(app);

    braid_Destroy(core);
    MPI_Finalize();

    return (0);
}
