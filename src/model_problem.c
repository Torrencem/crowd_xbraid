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
 * Example:       trischur-adv-diff.c
 *
 * Interface:     C
 *
 * Requires:      only C-language support
 *
 * Compile with:  make trischur-adv-diff
 *
 * Description:   Solves a linear optimal control problem in time-parallel:
 *
 *                min   0.5\int_0^T \int_0^1 (u(x,t)-u0(x))^2+alpha v(x,t)^2
 *dxdt
 *
 *                s.t.  du/dt + du/dx - nu d^2u/dx^2 = v(x,t)
 *                      u(0,t)=u(1,t)=0
 *                      u(x,0)=u0(x)
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

#include "lapacke.h"
#include "utils.c"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct {
    int myid;     /* Rank of the processor */
    double alpha; /* Relaxation parameter for objective function, v(x,t) */
    double nu;    /* Diffusion coefficent, which we take to be large */
    int ntime;    /* Total number of time-steps (starting at time 0) */
    int mspace;   /* Number of space points included in our state vector */
                  /* So including boundaries we have M+2 space points */
    double dx;    /* Spatial mesh spacing */

    int ilower;  /* Lower index for my proc */
    int iupper;  /* Upper index for my proc */
    int npoints; /* Number of time points on my proc */

    double **w;  /* Adjoint vectors at each time point on my proc */
    double *u0;  /* Initial condition */
    double *scr; /* Scratch space (enough for 3 spatial vectors) */
} my_App;

/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct {
    Vector u;
    Vector v;
    Vector w;
} my_Vector;

double left_boundary_solution_u(const double t) {
    return 0.0;
}

double right_boundary_solution_u(const double t) {
    return 1.0;
}

double left_boundary_solution_v(const double t) {
    return 0.0;
}

double right_boundary_solution_v(const double t) {
    return 0.0;
}

void compute_W_matrix(const double *w, const int n, Vector *WL, Vector *W, Vector *WU) {
    // Set WL
    *WL = zero_vector(n - 1);
    for (int i = 0; i < n - 1; i++) {
        (*WL)[i] = -w[i + 1];
    }

    // Set W
    *W = zero_vector(n);
    // Middle indices
    for (int i = 1; i < n - 2; i++) {
        (*W)[i] = 2 + w[i - 1] - w[i + 1];
    }
    // Edge cases
    (*W)[0] = -w[1] + 2;
    // Last element of W
    (*W)[n - 1] = 2 + w[n - 2];

    // Set WU
    *WU = zero_vector(n - 1);
    for (int i = 0; i < n - 1; i++) {
        (*WL)[i] = w[i];
    }
}

Vector compute_b_vector(const int n, const double *w, const double t) {
    double vleft = left_boundary_solution_v(t);
    double vright = right_boundary_solution_v(t);
    double wleft = w[0];
    double wright = w[n - 1];
    Vector ret = zero_vector(n);
    ret[0] = vleft * wleft;
    ret[n - 1] = - vright * wright;
    return ret;
}

// u must be allocated before calling
Vector apply_Phi(my_App *app, my_Vector *vec, double t, double dt) {
    Vector u_new = zero_vector(app->npoints);
    double beta = (dt / (2.0 * app->dx));
    // Space Boundary Conditions
    // j = 0
    u_new[0] = vec->u[0] +
        beta *
        ((vec->u[0]) * (vec->v[1] - left_boundary_solution_v(t)) +
         (vec->v[0]) * (vec->u[1] - left_boundary_solution_u(t)));
    // j = app->npoints
    int n = app->npoints - 1;
    u_new[n] = vec->u[n] +
        beta *
        ((vec->u[n]) * (right_boundary_solution_v(t) - vec->v[n - 1]) +
         (vec->v[n]) * (right_boundary_solution_u(t) - vec->u[n - 1]));
    // Non-Boundary points
    for (int j = 1; j < n; j++) {
        double uprev = vec->u[j];
        u_new[j] = uprev +
            beta * 
            ((uprev) * (vec->v[j + 1] - vec->v[j - 1]) +
             (vec->v[j]) * (vec->u[j + 1] - vec->u[j - 1]));
    }
    return u_new;
}


/*------------------------------------*/

/* Compute A(u) - f */

int apply_TriResidual(my_App *app, my_Vector *uleft, my_Vector *uright,
                      my_Vector *f, my_Vector *r, double dt) {
    double nu = (app->nu);
    double alpha = (app->alpha);
    int mspace = (app->mspace);
    double dx = (app->dx);
    double *u0 = (app->u0);
    double *scr = (app->scr);

    double *rtmp, *utmp;

    /* Create temporary vectors */
    rtmp = zero_vector(mspace);
    utmp = zero_vector(mspace);

    /* Compute action of center block */

    /* rtmp = U_i^{-1} u */
    vec_copy(mspace, (r->values), utmp);
    apply_Uinv(dt, dx, mspace, utmp);
    vec_copy(mspace, utmp, rtmp);

    /* rtmp = rtmp + D_i V_i^{-1} D_i^T u */
    vec_copy(mspace, (r->values), utmp);
    apply_DAdjoint(dt, dx, nu, mspace, utmp, scr);
    apply_Vinv(dt, dx, alpha, mspace, utmp);
    apply_D(dt, dx, nu, mspace, utmp, scr);
    vec_axpy(mspace, 1.0, utmp, 1.0, rtmp);

    /* rtmp = rtmp + Phi_i U_{i-1}^{-1} Phi_i^T u */
    /* This term is zero at time 0, since Phi_0 = 0 */
    if (uleft != NULL) {
        vec_copy(mspace, (r->values), utmp);
        apply_PhiAdjoint(dt, dx, nu, mspace, utmp, scr);
        apply_Uinv(dt, dx, mspace, utmp);
        apply_Phi(dt, dx, nu, mspace, utmp, scr);
        vec_axpy(mspace, 1.0, utmp, 1.0, rtmp);
    }

    /* Compute action of west block */
    if (uleft != NULL) {
        vec_copy(mspace, (uleft->values), utmp);
        apply_Uinv(dt, dx, mspace, utmp);
        apply_Phi(dt, dx, nu, mspace, utmp, scr);
        vec_axpy(mspace, -1.0, utmp, 1.0, rtmp);
    }

    /* Compute action of east block */
    if (uright != NULL) {
        vec_copy(mspace, (uright->values), utmp);
        apply_PhiAdjoint(dt, dx, nu, mspace, utmp, scr);
        apply_Uinv(dt, dx, mspace, utmp);
        vec_axpy(mspace, -1.0, utmp, 1.0, rtmp);
    }

    /* No change for index 0 */
    vec_copy(mspace, u0, utmp);
    vec_axpy(mspace, -1.0, utmp, 1.0, rtmp);
    vec_copy(mspace, u0, utmp);
    apply_Phi(dt, dx, nu, mspace, utmp, scr);
    vec_axpy(mspace, 1.0, utmp, 1.0, rtmp);

    /* Subtract rhs f */
    if (f != NULL) {
        /* rtmp = rtmp - f */
        vec_axpy(mspace, -1.0, (f->values), 1.0, rtmp);
    }

    /* Copy temporary residual vector into residual */
    vec_copy(mspace, rtmp, (r->values));

    /* Destroy temporary vectors */
    vec_destroy(rtmp);
    vec_destroy(utmp);

    return 0;
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

    /* Compute residual */
    apply_TriResidual(app, uleft, uright, f, r, dt);

    return 0;
}

/*------------------------------------*/

/* Solve A(u) = f */

int my_TriSolve(braid_App app, braid_Vector uleft, braid_Vector uright,
                braid_Vector fleft, braid_Vector fright, braid_Vector f,
                braid_Vector u, braid_TriStatus status) {
    double alpha = (app->alpha);
    double dx = (app->dx);

    double t, tprev, tnext, dt;
    double *utmp, scale;
    int iter;

    /* Get the time-step size */
    braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }

    // Compute u
    Vector u_new = apply_Phi(app, uleft, t, dt);

    // Compute v
    Vector RHS = zero_vector(app->mspace);
    vec_axpy(app->mspace, -1.0, uright->w, 1.0, RHS);
    vec_axpy(app->mspace, 2*dx*dt, u_new, 1.0, RHS);
    vec_axpy(app->mspace, -1.0, compute_b_vector(app->mspace, uright->w, t), 1.0, RHS);
    // Matrix

    // Compute w
    Vector w_new = zero_vector(app->mspace);
    vec_axpy(app->mspace, -2*dx*dt, v_new, 1.0, w_new);
    //Matrix

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
    u->values = zero_vector(mspace);

    for (i = 0; i <= mspace - 1; i++) {
        u->values[i] = ((double)braid_Rand()) / braid_RAND_MAX;
        //      u->values[i] = 1.0;
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
    u->values = zero_vector(mspace);
    vec_copy(mspace, (u->values), (v->values));

    *v_ptr = v;

    return 0;
}

/*------------------------------------*/

int my_Free(braid_App app, braid_Vector u) {
    free(u->values);
    free(u);

    return 0;
}

/*------------------------------------*/

int my_Sum(braid_App app, double alpha, braid_Vector x, double beta,
           braid_Vector y) {
    vec_axpy((app->mspace), alpha, (x->values), beta, (y->values));

    return 0;
}

/*------------------------------------*/

int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr) {
    int i;
    double dot = 0.0;
    int mspace = (app->mspace);
    for (i = 0; i < mspace; i++) {
        dot += (u->values)[i] * (u->values)[i];
    }
    *norm_ptr = sqrt(dot);

    return 0;
}

/*------------------------------------*/

// ZTODO: Need to compute u from adjoint and it reqires communication

int my_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus) {
    int done, index, ii;
    int mspace = (app->mspace);

    /* Print solution to file if simulation is over */
    braid_AccessStatusGetDone(astatus, &done);
    if (done) {
        braid_AccessStatusGetILowerUpper(astatus, &(app->ilower),
                                         &(app->iupper));
        (app->npoints) = (app->iupper) - (app->ilower) + 1;

        /* Allocate w array in app */
        if ((app->w) == NULL) {
            (app->w) = (double **)calloc((app->npoints), sizeof(double *));
        }

        braid_AccessStatusGetTIndex(astatus, &index);
        ii = index - (app->ilower);
        if (app->w[ii] != NULL) {
            free(app->w[ii]);
        }
        app->w[ii] = zero_vector(mspace);
        vec_copy(mspace, (u->values), (app->w[ii]));
    }

#if 0
   //  Below prints U, V, and W after selected iterations. This can then be
   //  plotted to show how the space-time solution changes after iterations.
   char  filename[255];
   FILE *file;
   int  iter;
   braid_AccessStatusGetIter(astatus, &iter);
   braid_AccessStatusGetTIndex(astatus, &index);
   /* file format is advec-diff-btcs.out.{iteration #}.{time index} */
   if(iter%1==0){
      sprintf(filename, "%s.%04d.%04d", "out/advec-diff-btcs.v.out", iter, index);
      file = fopen(filename, "w");
      for(i = 0; i<mspace; i++){
         if(i<mspace-1){
            fprintf(file, "%1.14e, ", (u->values)[i]);
         }
         else{
            fprintf(file, "%1.14e", (u->values)[i]);
         }
      }
      fflush(file);
      fclose(file);
   }
#endif

    return 0;
}

/*------------------------------------*/

int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus bstatus) {
    *size_ptr = (app->mspace) * sizeof(double);
    return 0;
}

/*------------------------------------*/

int my_BufPack(braid_App app, braid_Vector u, void *buffer,
               braid_BufferStatus bstatus) {
    double *dbuffer = buffer;

    vec_copy((app->mspace), (u->values), dbuffer);
    braid_BufferStatusSetSize(bstatus, (app->mspace) * sizeof(double));

    return 0;
}

/*------------------------------------*/

int my_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr,
                 braid_BufferStatus bstatus) {
    my_Vector *u = NULL;
    double *dbuffer = buffer;

    /* Allocate memory */
    u = (my_Vector *)malloc(sizeof(my_Vector));
    u->values = zero_vector((app->mspace));

    /* Unpack the buffer */
    vec_copy((app->mspace), dbuffer, (u->values));

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
    double alpha, nu;
    int max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
    int access_level, print_level;
    double tol;
    double time;
    double *u0, *scr;
    int i;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Define space domain. Space domain is between 0 and 1, mspace defines the
     * number of steps */
    mspace = 8;
    ntime = 256;

    /* Define some optimization parameters */
    alpha = .005; /* parameter in the objective function */
    nu = .03;     /* parameter in PDE (used 0.3 in RIPS) */

    /* Define some Braid parameters */
    max_levels = 30;
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
            printf("  -nu <nu>                : Constant Parameter in PDE  \n");
            printf("  -alpha <alpha>          : Constant Parameter in "
                   "Objective Function  \n");
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
        } else if (strcmp(argv[arg_index], "-nu") == 0) {
            arg_index++;
            nu = atof(argv[arg_index++]);
        } else if (strcmp(argv[arg_index], "-alpha") == 0) {
            arg_index++;
            alpha = atof(argv[arg_index++]);
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

    /* Set up initial condition */
    u0 = zero_vector(mspace);
    for (i = 0; i < mspace / 2; i++) {
        u0[i] = 1;
    }
    for (i = mspace / 2; i < mspace; i++) {
        u0[i] = 0;
    }

    /* Set up scratch space */
    scr = zero_vector(3 * mspace);

    /* Set up the app structure */
    app = (my_App *)malloc(sizeof(my_App));
    app->myid = rank;
    app->ntime = ntime;
    app->mspace = mspace;
    app->dx = dx;
    app->nu = nu;
    app->alpha = alpha;
    app->w = NULL;
    app->u0 = u0;
    app->scr = scr;

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

    /* Writes solutions to files */
    if (access_level > 0) {
        char filename[255];
        FILE *file;
        int i, j;
        double **w = (app->w);
        double **v;
        double **u;

        /* Print adjoint w to file */
        sprintf(filename, "%s.%03d", "trischur-adv-diff.out.w", (app->myid));
        file = fopen(filename, "w");
        for (i = 0; i < (app->npoints); i++) {
            fprintf(file, "%05d: ", ((app->ilower) + i + 1));
            for (j = 0; j < mspace; j++) {
                if (j == mspace - 1) {
                    fprintf(file, "% 1.14e", w[i][j]);
                } else {
                    fprintf(file, "% 1.14e, ", w[i][j]);
                }
            }
            fprintf(file, "\n");
        }
        fflush(file);
        fclose(file);

        /* Compute control v from adjoint w and print to file */
        sprintf(filename, "%s.%03d", "trischur-adv-diff.out.v", (app->myid));
        file = fopen(filename, "w");
        v = (double **)malloc(app->npoints * sizeof(double *));
        for (i = 0; i < app->npoints; i++) {
            v[i] = zero_vector((app->mspace));
        }
        for (i = 0; i < (app->npoints); i++) {
            vec_copy(mspace, w[i], v[i]);
            apply_DAdjoint(dt, dx, nu, mspace, v[i], scr);
            apply_Vinv(dt, dx, alpha, mspace, v[i]);

            fprintf(file, "%05d: ", ((app->ilower) + i + 1));
            for (j = 0; j < (app->mspace); j++) {
                if (j == mspace - 1) {
                    fprintf(file, "% 1.14e", v[i][j]);
                } else {
                    fprintf(file, "% 1.14e, ", v[i][j]);
                }
            }
            fprintf(file, "\n");
        }
        fflush(file);
        fclose(file);

        /* Compute state u from adjoint w and print to file */
        sprintf(filename, "%s.%03d", "trischur-adv-diff.out.u", (app->myid));
        file = fopen(filename, "w");
        u = (double **)malloc(app->npoints * sizeof(double *));
        for (i = 0; i < app->npoints; i++) {
            u[i] = zero_vector((app->mspace));
        }
        for (i = 0; i < (app->npoints); i++) {
            if (i != app->npoints - 1) {
                vec_copy(mspace, w[i + 1], u[i]);
                apply_PhiAdjoint(dt, dx, nu, mspace, u[i], scr);
                vec_axpy(mspace, -1.0, w[i], 1.0, u[i]);
                vec_axpy(mspace, dx * dt, u0, 1.0, u[i]);
                apply_Uinv(dt, dx, mspace, u[i]);
            } else {
                vec_copy(mspace, w[i], u[i]);
                vec_scale(mspace, -1.0, u[i]);
                vec_axpy(mspace, dx * dt, u0, 1.0, u[i]);
                apply_Uinv(dt, dx, mspace, u[i]);
            }

            fprintf(file, "%05d: ", ((app->ilower) + i + 1));
            for (j = 0; j < mspace; j++) {
                if (j == mspace - 1) {
                    fprintf(file, "% 1.14e", u[i][j]);
                } else {
                    fprintf(file, "% 1.14e, ", u[i][j]);
                }
            }
            fprintf(file, "\n");
        }
        fflush(file);
        fclose(file);

        /* Print initial guess u0 to file */
        sprintf(filename, "%s.%03d", "trischur-adv-diff.out.u0", (app->myid));
        file = fopen(filename, "w");
        for (j = 0; j < mspace; j++) {
            if (j != mspace - 1) {
                fprintf(file, "% 1.14e, ", u0[j]);
            } else {
                fprintf(file, "% 1.14e", u0[j]);
            }
        }
        fflush(file);
        fclose(file);

#if 0
      // Calculates value of objective function
      double objective_val=0;

      for(i=0; i<ntime; i++)
      {
         for(int j=0; j<mspace; j++)
         {
            objective_val += ((u[i][j] - u0[j]) * (u[i][j] - u0[j]) +
                              alpha*v[i][j]*v[i][j] ) * dx;
         }
         objective_val *= dt;
      }
      printf("Objective Function Value: %f \n", objective_val);
#endif

        for (i = 0; i < app->npoints; i++) {
            vec_destroy(w[i]);
            vec_destroy(v[i]);
            vec_destroy(u[i]);
        }
        free(w);
        free(v);
        free(u);
    }

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

    vec_destroy(app->u0);
    vec_destroy(app->scr);
    free(app);

    braid_Destroy(core);
    MPI_Finalize();

    return (0);
}
