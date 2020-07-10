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
 * Example:       model_problem.cpp
 *
 * Interface:     C
 *
 * Requires:      Lapacke
 *
 * Compile with:  make model_problem
 *
 * Description:   Solves a multilinear optimal control problem in time-parallel:
 *
 *                min   \int_0^T \int_0^1 u(x, t)^2+v(x, t)^2 dx dt
 *dxdt
 *
 *                s.t.  du/dt = c * d/dx(uv) where c is a constant
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "braid.h"
#include "braid_test.h"
#define PI 3.14159265
#define g(dt, dx) dt / (2 * dx)
//#define g(dt,dx) 0.0
#define b(dt, dx, nu) nu *dt / (dx * dx)
//#define c 1.0

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

    int ilower;  /* Lower index for my proc */
    int iupper;  /* Upper index for my proc */
    int npoints; /* Number of time points on my proc */
    FILE *ufile, *vfile, *wfile;
} my_App;

/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct {
    Vector u;
    Vector v;
    Vector w;
} my_Vector;

/*
double sign (const double c){ // Time steps and space steps are 
    if(c < 0.0){
        // If c is negative, for the negative beta we want the minimum of (c, 0)
        return (c > 0.0) ? c : 0.0;
    }else if(c > 0.0){
        // Else if c is positive we want the maximum of (c, 0)
        return (c > 0.0) ? 0.0 : c;
    }
}
*/

double left_boundary_solution_u(const double t) { return 0.0; }

double right_boundary_solution_u(const double t) { return 1.0; }

double left_boundary_solution_v(const double t) { return 0.0; }

double right_boundary_solution_v(const double t) { return 0.0; }

double initial_condition_u(const double x) { return 0.25 * x; }

double initial_condition_v(const double x) { return 0.0; }

void compute_L_matrix(const my_App *app, const double *v, const int n,
                      Matrix *L, const double dt, const double t) {
    Vector LL, LC, LU;
    // Set LL
    double c = 1.0;
    double beta_neg = fmin(c, 0.0) * (dt / (app->dx));
    double beta_pos = fmax(c, 0.0) * (dt / (app->dx));

    LL = zero_vector(n - 1);
    for (int i = 0; i < n - 1; i++) {
        LL[i] = -beta_pos * v[i];
    }

    // Set L
    LC = zero_vector(n);

    for(int i = 0; i < n; i++){
        LC[i] = 1 + beta_pos * (v[i - 1] - 2*v[i]) - 
                       beta_neg * (2*v[i] - v[i + 1]);
    }

    // Set LU
    LU = zero_vector(n - 1);

    for (int i = 0; i < n - 1; i++){
        LU[0] = beta_neg * v[i];
    }

    *L = tridiag_to_matrix(LL, LC, LU, n);

    set_element(*L, n, n, 0, n - 1, -beta_pos * v[0]);
    set_element(*L, n, n, n - 1, 0, beta_neg * v[n - 1]);
}

Vector compute_b_vector(const int n, const double *w, const double t) {
    double vleft = left_boundary_solution_v(t);
    double vright = right_boundary_solution_v(t);
    double wleft = w[0];
    double wright = w[n - 1];
    Vector ret = zero_vector(n);
    ret[0] = vleft * wleft;
    ret[n - 1] = -vright * wright;
    return ret;
}

Vector apply_Phi(const my_App *app, const Vector u, const Vector v,
                 const double t, const double dt) {
    Vector u_new = zero_vector(app->mspace);
    double c = 1.0;
    double beta_neg = fmin(c,0.0) * (dt / (app->dx));
    double beta_pos = fmax(c,0.0) * (dt / (app->dx));
    int n = app->mspace-1;
    // Space Boundary Conditions
    // j = 0
    u_new[0] = u[0] * (1 - beta_neg * (2*v[0] - v[1]) + //is right boundary solution v == v[n] and can we do that in this?
                              beta_pos * (v[n - 1] - 2*v[0])) + beta_neg * u[1] * v[0] - beta_pos * v[0] * u[n - 1]; //periodic boundaries
   
    // j = app->mspace
    u_new[n - 1] =
        u[n - 1] * (1 - beta_neg * (2 * v[n - 1] - v[0]) + beta_pos * (v[n - 2] - 2 * v[n - 1])) - beta_pos * v[n - 1] * u[n - 2] + beta_neg * v[n - 1] * u[0];

    // Non-Boundary points
    for (int j = 1; j < n - 1; j++) {
        double uprev = u[j];
        u_new[j] = uprev - beta_neg * ((uprev) * (v[j + 1] - v[j]) + v[j] * (u[j + 1] - uprev)) - 
                               beta_pos * (uprev * (v[j] - v[j - 1]) + v[j] * (uprev - u[j - 1]));
    }
    return u_new;
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

    Vector u_res;
    // Compute u
    if (uleft == NULL) {
        // Collect initial conditions
        Vector u = zero_vector(app->mspace);
        Vector v = zero_vector(app->mspace);
        for (int x = 0; x < app->mspace; x++) {
            u[x] = initial_condition_u(x);
            v[x] = initial_condition_v(x);
        }
        u_res = apply_Phi(app, u, v, t, dt);
    } else {
        u_res = apply_Phi(app, uleft->u, uleft->v, t, dt);
    }
    vec_axpy(app->mspace, 1.0, r->u, -1.0, u_res);

    Vector v_res = zero_vector(app->mspace);
    if (uright != NULL) {
        vec_copy(app->mspace, uright->w, v_res);
        Vector L = NULL;
        compute_L_matrix(app, r->u, app->mspace, &L, dt, t);
        matmul(L, app->mspace, app->mspace, &v_res);
    }
    vec_axpy(app->mspace, 2.0 * app->dx * dt, r->v, -1.0, v_res);

    Vector w_res = zero_vector(app->mspace);
    if (uright != NULL) {
        vec_copy(app->mspace, uright->w, w_res);
        Vector L = NULL;
        compute_L_matrix(app, r->v, app->mspace, &L, dt, t);
        for (int i = 0; i < app->mspace; i++) {
            double val = get_element(L, app->mspace, app->mspace, i, i);
            set_element(L, app->mspace, app->mspace, i, i, val);
        }
        matmul(L, app->mspace, app->mspace, &w_res);
    }
    vec_axpy(app->mspace, -2.0 * app->dx * dt, r->u, 1.0, w_res);

    vec_copy(app->mspace, u_res, r->u);
    vec_copy(app->mspace, v_res, r->v);
    vec_copy(app->mspace, w_res, r->w);

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

    Vector u_new;
    // Compute u
    if (uleft == NULL) {
        // Collect initial conditions
        Vector u = zero_vector(app->mspace);
        Vector v = zero_vector(app->mspace);
        for (int x = 0; x < app->mspace; x++) {
            u[x] = initial_condition_u(x);
            v[x] = initial_condition_v(x);
        }
        u_new = apply_Phi(app, u, v, t, dt);
    } else {
        u_new = apply_Phi(app, uleft->u, uleft->v, t, dt);
    }

    // Compute v
    Vector v_new = zero_vector(app->mspace);
    if (uright != NULL) {
        Vector L = NULL;
        compute_L_matrix(app, u_new, app->mspace, &L, dt, t);
        vec_copy(app->mspace, uright->w, v_new);
        matmul(L, app->mspace, app->mspace, &v_new);
        vec_scale(app->mspace, 1.0 / (2.0 * dx * dt), v_new);
    }
    // Compute w
    Vector w_new = zero_vector(app->mspace);
    if (uright != NULL) {
        vec_copy(app->mspace, uright->w, w_new);
        // Compute I+Lv
        Vector L = NULL;
        compute_L_matrix(app, v_new, app->mspace, &L, dt, t);
        for (int i = 0; i < app->mspace; i++) {
            double val = get_element(L, app->mspace, app->mspace, i, i);
            set_element(L, app->mspace, app->mspace, i, i, val + 1.0);
        }
        matmul(L, app->mspace, app->mspace, &w_new);
    }
    vec_axpy(app->mspace, -2.0 * dx * dt, u_new, 1.0, w_new);

    vec_copy(app->mspace, u_new, u->u);
    vec_copy(app->mspace, v_new, u->v);
    vec_copy(app->mspace, w_new, u->w);

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
    u->u = zero_vector(mspace);
    u->v = zero_vector(mspace);
    u->w = zero_vector(mspace);

    for (i = 0; i <= mspace - 1; i++) {
        u->u[i] = ((double)braid_Rand()) / braid_RAND_MAX;
        u->v[i] = ((double)braid_Rand()) / braid_RAND_MAX;
        u->w[i] = ((double)braid_Rand()) / braid_RAND_MAX;
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

    v->u = zero_vector(app->mspace);
    v->v = zero_vector(app->mspace);
    v->w = zero_vector(app->mspace);

    vec_copy(mspace, (u->u), (v->u));
    vec_copy(mspace, (u->v), (v->v));
    vec_copy(mspace, (u->w), (v->w));
    *v_ptr = v;

    return 0;
}

/*------------------------------------*/

int my_Free(braid_App app, braid_Vector u) {
    free(u->u);
    free(u->v);
    free(u->w);
    free(u);

    return 0;
}

/*------------------------------------*/

int my_Sum(braid_App app, double alpha, braid_Vector x, double beta,
           braid_Vector y) {

    vec_axpy((app->mspace), alpha, (x->u), beta, (y->u));
    vec_axpy((app->mspace), alpha, (x->v), beta, (y->v));
    vec_axpy((app->mspace), alpha, (x->w), beta, (y->w));
    return 0;
}

/*------------------------------------*/

int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr) {
    int i;
    double dot = 0.0;
    int mspace = (app->mspace);
    for (i = 0; i < mspace; i++) {
        dot += (u->u)[i] * (u->u)[i];
        dot += (u->v)[i] * (u->v)[i];
        dot += (u->w)[i] * (u->w)[i];
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

        fprintf((app->ufile), "%05d: ", index);
        fprintf((app->vfile), "%05d: ", index);
        fprintf((app->wfile), "%05d: ", index);
        for (j = 0; j < (mspace - 1); j++) {
            fprintf((app->ufile), "% 1.14e, ", (u->u[j]));
            fprintf((app->vfile), "% 1.14e, ", (u->v[j]));
            fprintf((app->wfile), "% 1.14e, ", (u->w[j]));
        }
        fprintf((app->ufile), "% 1.14e\n", (u->u[j]));
        fprintf((app->vfile), "% 1.14e\n", (u->v[j]));
        fprintf((app->wfile), "% 1.14e\n", (u->w[j]));
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

    vec_copy((app->mspace), (u->u), dbuffer);
    dbuffer += (app->mspace);
    vec_copy((app->mspace), (u->v), dbuffer);
    dbuffer += (app->mspace);
    vec_copy((app->mspace), (u->w), dbuffer);
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
    vec_create((app->mspace), &(u->u));
    vec_create((app->mspace), &(u->v));
    vec_create((app->mspace), &(u->w));

    /* Unpack the buffer */
    vec_copy((app->mspace), dbuffer, (u->u));
    dbuffer += (app->mspace);
    vec_copy((app->mspace), dbuffer, (u->v));
    dbuffer += (app->mspace);
    vec_copy((app->mspace), dbuffer, (u->w));

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
    sprintf(filename, "%s.%03d", "modelproblem.out.u", (app->myid));
    (app->ufile) = fopen(filename, "w");
    sprintf(filename, "%s.%03d", "modelproblem.out.v", (app->myid));
    (app->vfile) = fopen(filename, "w");
    sprintf(filename, "%s.%03d", "modelproblem.out.w", (app->myid));
    (app->wfile) = fopen(filename, "w");

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

    //   {
    //      char    filename[255];
    //      FILE   *file;
    //
    //      //Note that this out file appends the number of time steps
    //      sprintf(filename, "%s.%d", "trischur-adv-diff.time", ntime);
    //
    //      file = fopen(filename, "w");
    //      fprintf(file, "%f", time);
    //      fflush(file);
    //      fclose(file);
    //   }

    //   /* RDF Testing */
    //   {
    //      double     cscale = (1/(dx*dt))*(2 + dt*dt/alpha);
    //      double     lscale = (1/(dx*dt))*(1 + dt*dt/alpha);
    //      my_Vector *e = (my_Vector *) malloc(sizeof(my_Vector));
    //      my_Vector *z = (my_Vector *) malloc(sizeof(my_Vector));
    //
    //      (e->values) = (double*) malloc( mspace*sizeof(double) );
    //      (z->values) = (double*) malloc( mspace*sizeof(double) );
    //      for(i = 0; i < mspace; i++)
    //      {
    //         (e->values[i]) = 1.0;
    //         (z->values[i]) = 0.0;
    //      }
    //      apply_TriResidual(app, z, z, NULL, e, 1, dt);
    //
    //      for(i = 0; i < mspace; i++)
    //      {
    //         (e->values[i]) = 1.0;
    //         (z->values[i]) = 0.0;
    //      }
    //      apply_TriResidual(app, NULL, z, NULL, e, 1, dt);
    //
    //      for(i = 0; i < mspace; i++)
    //      {
    //         (e->values[i]) = 1.0;
    //         (z->values[i]) = 0.0;
    //      }
    //      apply_TriResidual(app, z, NULL, NULL, e, 1, dt);
    //   }

    free(app);

    braid_Destroy(core);
    MPI_Finalize();

    return (0);
}