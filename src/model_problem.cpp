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
 * Description:   Solves a multilinear optimal control problem in time-parallel:
 *
 *                min   \int_0^T \int_0^1 u(x, t)^2+v(x, t)^2 dx dt
 *dxdt
 *
 *                s.t.  du/dt = d/dx(uv)
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "tribraid.hpp"
#include "braid_test.h"
#define PI 3.14159265
#define g(dt, dx) dt / (2 * dx)
//#define g(dt,dx) 0.0
#define b(dt, dx, nu) nu *dt / (dx * dx)

#include "lapacke.h"
#include "utils.cpp"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

class BraidVector
{
    public:
    Vector u, v, w;

    BraidVector(Vector u_, Vector v_, Vector w_) : u(u_), v(v_), w(w_) { }

    virtual ~BraidVector() { };
};

double left_boundary_solution_u(const double t) { return 0.0; }

double right_boundary_solution_u(const double t) { return 1.0; }

double left_boundary_solution_v(const double t) { return 0.0; }

double right_boundary_solution_v(const double t) { return 0.0; }

double initial_condition_u(const double x) { return 0.25 * x; }

double initial_condition_v(const double x) { return 0.0; }

class MyBraidApp : public TriBraidApp {
protected:
    // BraidApp defines tstart, tstop, ntime and comm_t
public:
    int ntime, mspace, ilower, iupper, npoints, myid;
    double dx;
    FILE *ufile, *vfile, *wfile;

    // Main constructor
    MyBraidApp(MPI_Comm comm_t__, int rank_, double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);

    virtual ~MyBraidApp() {};

    // Define the Vraid Wrapper routines
    // Note: braid_Vector is the type BraidVector*
    virtual int Clone(braid_Vector u_,
                      braid_Vector *v_ptr) override;

    virtual int Init(double t,
                     braid_Vector *v_ptr) override;

    virtual int Free(braid_Vector u_) override;

    virtual int Sum(double alpha,
                    braid_Vector x_,
                    double beta,
                    braid_Vector y_) override;

    virtual int SpatialNorm(braid_Vector u_,
                            double *norm_ptr) override;

    virtual int Access(braid_Vector       u_,
                            BraidAccessStatus &astatus) override;

    virtual int TriResidual(braid_Vector uleft_,
                            braid_Vector uright_,
                            braid_Vector f_,
                            braid_Vector r_,
                            BraidTriStatus &pstatus) override;

    virtual int TriSolve(braid_Vector uleft_,
                         braid_Vector uright_,
                         braid_Vector fleft_,
                         braid_Vector fright_,
                         braid_Vector f_,
                         braid_Vector u_,
                         BraidTriStatus &pstatus) override;

    virtual int BufSize(braid_Int         *size_ptr,
                             BraidBufferStatus &bstatus) override;

    virtual int BufPack(braid_Vector       u_,
                             void              *buffer,
                             BraidBufferStatus &bstatus) override;

    virtual int BufUnpack(void              *buffer,
                               braid_Vector      *u_ptr,
                               BraidBufferStatus &bstatus) override;

    // Not needed
    virtual int Residual(braid_Vector u_,
                            braid_Vector r_,
                            BraidStepStatus &pstatus) override { return 0; }

    const Tridiag_Matrix compute_L_matrix(Vector &V, const double dt, const double t);

    const Vector apply_Phi(Vector &u, Vector &v,
                 const double t, const double dt);
};

MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_, double tstop_, int ntime_) : TriBraidApp(comm_t_, tstart_, tstop_, ntime_) {

}

const Tridiag_Matrix MyBraidApp::compute_L_matrix(Vector &v, const double dt, const double t) {
    int n = v->len;
    Tridiag_Matrix L = new Tridiag_Matrix_(n);
    for (int i = 0; i < n - 1; i++) {
        L->al[i] = (*v)[i];
    }

    for (int i = 1; i < n - 2; i++) {
#ifdef MODEL_BACKWARDS
        L->a[i] = (*v)[i] - (*v)[i - 1];
#else
        L->a[i] = (*v)[i + 1] - (*v)[i - 1];
#endif
    }

    L->a[0] = (*v)[1] - left_boundary_solution_v(t);

    L->a[n - 1] = right_boundary_solution_v(t) - (*v)[n - 2];

    for (int i = 0; i < n - 1; i++) {
        L->al[i] = -(*v)[i];
    }

    return L;
}

Vector compute_b_vector(Vector &w, const double t) {
    int n = w->len;
    double vleft = left_boundary_solution_v(t);
    double vright = right_boundary_solution_v(t);
    double wleft = (*w)[0];
    double wright = (*w)[n - 1];
    Vector ret = new Vector_(n);
    ret[0] = vleft * wleft;
    ret[n - 1] = -vright * wright;
    return ret;
}

const Vector MyBraidApp::apply_Phi(Vector &u_, Vector &v_,
                 const double t, const double dt) {
    Vector u_new = new Vector_(mspace);
#ifdef MODEL_BACKWARDS
    double beta = (dt / dx);
#else
    double beta = (dt / (2.0 * dx));
#endif

    Vector_ &u = *u_;
    Vector_ &v = *v_;

    // Space Boundary Conditions
    // j = 0
#ifdef MODEL_BACKWARDS
    (*u_new)[0] = u[0] + beta * ((u[0]) * (v[0] - left_boundary_solution_v(t)) +
                              (v[0]) * (u[0] - left_boundary_solution_u(t)));
#else
    (*u_new)[0] = u[0] + beta * ((u[0]) * (v[1] - left_boundary_solution_v(t)) +
                              (v[0]) * (u[1] - left_boundary_solution_u(t)));
#endif

    // j = app->mspace
#ifdef MODEL_BACKWARDS
    int n = mspace - 1;
    (*u_new)[n] =
        u[n] + beta * ((u[n]) * (v[n] - v[n - 1]) + (v[n]) * (u[n] - u[n - 1]));
#else
    int n = mspace - 1;
    (*u_new)[n] =
        u[n] + beta * ((u[n]) * (right_boundary_solution_v(t) - v[n - 1]) +
                       (v[n]) * (right_boundary_solution_u(t) - u[n - 1]));
#endif

    // Non-Boundary points
    for (int j = 1; j < n; j++) {
        double uprev = u[j];
#ifdef MODEL_BACKWARDS
        (*u_new)[j] = uprev + beta * ((uprev) * (v[j + 1] - v[j - 1]) +
                                   (v[j]) * (u[j + 1] - u[j - 1]));
#else
        (*u_new)[j] = uprev + beta * ((uprev) * (v[j] - v[j - 1]) +
                                   (v[j]) * (u[j] - u[j - 1]));
#endif
    }
    return u_new;
}

/*------------------------------------*/

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */

int MyBraidApp::TriResidual(braid_Vector uleft_, braid_Vector uright_,
                   braid_Vector f_, braid_Vector r_, BraidTriStatus &status) {
    BraidVector *uleft = (BraidVector *) uleft_;
    BraidVector *uright = (BraidVector *) uright_;
    // BraidVector *f = (BraidVector *) f_;
    BraidVector *r = (BraidVector *) r_;

    double t, tprev, tnext, dt;
    int level, index;

    status.GetTriT(&t, &tprev, &tnext);
    status.GetLevel(&level);
    status.GetTIndex(&index);

    /* Get the time-step size */
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }

    Vector u_res = new Vector_(1);
    // Compute u
    if (uleft == NULL) {
        // Collect initial conditions
        Vector u = new Vector_(mspace);
        Vector v = new Vector_(mspace);
        for (int x = 0; x < mspace; x++) {
            (*u)[x] = initial_condition_u(x);
            (*v)[x] = initial_condition_v(x);
        }
        u_res = this->apply_Phi(u, v, t, dt);
    } else {
        u_res = this->apply_Phi(uleft->u, uleft->v, t, dt);
    }
    u_res = *r->u - u_res;

    Vector v_res = new Vector_(mspace);
    if (uright != NULL) {
        Tridiag_Matrix L = this->compute_L_matrix(r->u, dt, t);
        v_res = *L * uright->w;
    }
    v_res = *(*r->v * (2.0 * dx * dt)) - v_res;

    Vector w_res = new Vector_(mspace);
    if (uright != NULL) {
        auto L = compute_L_matrix(r->v, dt, t);
        (*L->a) += 1;
        w_res = *L * uright->w;
    }
    w_res = *(*r->u * (-2.0 * dx * dt)) + w_res;

    r->u = u_res;
    r->v = v_res;
    r->w = w_res;

    return 0;
}

/*------------------------------------*/

/* Solve A(u) = f */

int MyBraidApp::TriSolve(braid_Vector uleft_, braid_Vector uright_,
                braid_Vector fleft_, braid_Vector fright_, braid_Vector f_,
                braid_Vector u_, BraidTriStatus &status) {
    BraidVector *uleft = (BraidVector *) uleft_;
    BraidVector *uright = (BraidVector *) uright_;
    // BraidVector *f = (BraidVector *) f_;
    BraidVector *u = (BraidVector *) u_;
    // BraidVector *fleft = (BraidVector *) fleft_;
    // BraidVector *fright = (BraidVector *) fright_;

    double t, tprev, tnext, dt;

    /* Get the time-step size */
    status.GetTriT(&t, &tprev, &tnext);
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }

    Vector u_new = new Vector_(1);
    // Compute u
    if (uleft == NULL) {
        // Collect initial conditions
        Vector u = new Vector_(mspace);
        Vector v = new Vector_(mspace);
        for (int x = 0; x < mspace; x++) {
            (*u)[x] = initial_condition_u(x);
            (*v)[x] = initial_condition_v(x);
        }
        u_new = this->apply_Phi(u, v, t, dt);
    } else {
        u_new = this->apply_Phi(uleft->u, uleft->v, t, dt);
    }

    // Compute v
    Vector v_new = new Vector_(mspace);
    if (uright != NULL) {
        auto L = this->compute_L_matrix(u_new, dt, t);
        v_new = *(*L * uright->w) * (1.0 / (2.0 * dx * dt));
    }

    // Compute w
    Vector w_new = new Vector_(mspace);
    if (uright != NULL) {
        auto L = this->compute_L_matrix(v_new, dt, t);
        w_new = *L * uright->w;
    }
    w_new = *(*u_new * (-2.0 * dx * dt)) + w_new;

    u->u = u_new;
    u->v = u_new;
    u->w = u_new;

    /* no refinement */
    status.SetRFactor(1);

    return 0;
}

/*------------------------------------*/

/* This is only called from level 0 */

int MyBraidApp::Init(double t, braid_Vector *u_ptr_) {
    BraidVector **u_ptr = (BraidVector **) u_ptr_;
    int i;
    Vector uu = new Vector_(mspace);
    Vector uv = new Vector_(mspace);
    Vector uw = new Vector_(mspace);
    BraidVector *u = new BraidVector(uu, uv, uw);

    for (i = 0; i <= mspace - 1; i++) {
        (*u->u)[i] = ((double)braid_Rand()) / braid_RAND_MAX;
        (*u->v)[i] = ((double)braid_Rand()) / braid_RAND_MAX;
        (*u->w)[i] = ((double)braid_Rand()) / braid_RAND_MAX;
    }

    *u_ptr = u;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Clone(braid_Vector u_, braid_Vector *v_ptr_) {
    BraidVector *u = (BraidVector *) u_;
    BraidVector **v_ptr = (BraidVector **) v_ptr_;

    Vector uu = new Vector_(*u->u);
    Vector uv = new Vector_(*u->u);
    Vector uw = new Vector_(*u->u);

    BraidVector *v = new BraidVector(uu, uv, uw);
    *v_ptr = v;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Free(braid_Vector u_) {
    BraidVector *u = (BraidVector *) u_;
    delete u;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Sum(double alpha, braid_Vector x_, double beta,
           braid_Vector y_) {
    BraidVector *y = (BraidVector *) y_;
    BraidVector *x = (BraidVector *) x_;

    y->u = *(*x->u * alpha) + y->u;
    y->v = *(*x->v * alpha) + y->v;
    y->w = *(*x->w * alpha) + y->w;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::SpatialNorm(braid_Vector u_, double *norm_ptr) {
    BraidVector *u = (BraidVector *) u_;
    int i;
    double dot = 0.0;
    for (i = 0; i < mspace; i++) {
        dot += (*u->u)[i] * (*u->u)[i];
        dot += (*u->v)[i] * (*u->v)[i];
        dot += (*u->w)[i] * (*u->w)[i];
    }
    *norm_ptr = sqrt(dot);

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Access(braid_Vector u_, BraidAccessStatus &astatus) {
    BraidVector *u = (BraidVector *) u_;

    int done, index;

    /* Print solution to file if simulation is over */
    astatus.GetDone(&done);
    if (done) {
        int j;

        astatus.GetTIndex(&index);

        fprintf(ufile, "%05d: ", index);
        fprintf(vfile, "%05d: ", index);
        fprintf(wfile, "%05d: ", index);
        for (j = 0; j < (mspace - 1); j++) {
            fprintf(ufile, "% 1.14e, ", (*u->u)[j]);
            fprintf(vfile, "% 1.14e, ", (*u->v)[j]);
            fprintf(wfile, "% 1.14e, ", (*u->w)[j]);
        }
        fprintf(ufile, "% 1.14e\n", (*u->u)[j]);
        fprintf(vfile, "% 1.14e\n", (*u->v)[j]);
        fprintf(wfile, "% 1.14e\n", (*u->w)[j]);
    }

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufSize(int *size_ptr, BraidBufferStatus &bstatus) {
    *size_ptr = 3 * mspace * sizeof(double);
    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufPack(braid_Vector u_, void *buffer,
               BraidBufferStatus &bstatus) {
    BraidVector *u = (BraidVector *) u_;
    double *dbuffer = (double *) buffer;

    for (int i = 0; i < mspace; i++) {
        u->u[i] = dbuffer[i];
        dbuffer[i] = (*u->u)[i];
    }
    dbuffer += mspace;

    for (int i = 0; i < mspace; i++) {
        u->v[i] = dbuffer[i];
        dbuffer[i] = (*u->v)[i];
    }
    dbuffer += mspace;
    
    for (int i = 0; i < mspace; i++) {
        u->w[i] = dbuffer[i];
        dbuffer[i] = (*u->w)[i];
    }
    bstatus.SetSize(3 * mspace * sizeof(double));

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufUnpack(void *buffer, braid_Vector *u_ptr_,
                 BraidBufferStatus &bstatus) {
    BraidVector **u_ptr = (BraidVector **) u_ptr_;
    Vector uu = new Vector_(mspace);
    Vector uv = new Vector_(mspace);
    Vector uw = new Vector_(mspace);
    BraidVector *u = new BraidVector(uu, uv, uw);
    double *dbuffer = (double *) buffer;

    for (int i = 0; i < mspace; i++) {
        u->u[i] = dbuffer[i];
        dbuffer[i] = (*u->u)[i];
    }
    dbuffer += mspace;

    for (int i = 0; i < mspace; i++) {
        u->v[i] = dbuffer[i];
        dbuffer[i] = (*u->v)[i];
    }
    dbuffer += mspace;
    
    for (int i = 0; i < mspace; i++) {
        u->w[i] = dbuffer[i];
        dbuffer[i] = (*u->w)[i];
    }
    bstatus.SetSize(3 * mspace * sizeof(double));

    *u_ptr = u;
    return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    double tstart, tstop, dx, start, end;
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

    /* Set up the app structure */
    auto app = MyBraidApp(MPI_COMM_WORLD, rank, tstart, tstop, ntime);
    app.myid = rank;
    app.ntime = ntime;
    app.mspace = mspace;
    app.dx = dx;

    char filename[255];
    sprintf(filename, "%s.%03d", "modelproblem.out.u", (app.myid));
    (app.ufile) = fopen(filename, "w");
    sprintf(filename, "%s.%03d", "modelproblem.out.v", (app.myid));
    (app.vfile) = fopen(filename, "w");
    sprintf(filename, "%s.%03d", "modelproblem.out.w", (app.myid));
    (app.wfile) = fopen(filename, "w");

    /* Initialize XBraid */

    auto core = BraidTriCore(MPI_COMM_WORLD, &app);

    /* Set some XBraid(_Adjoint) parameters */
    core.SetMaxLevels(max_levels);
    core.SetMinCoarse(min_coarse);
    core.SetNRelax(-1, nrelax);
    if (max_levels > 1) {
        core.SetNRelax(max_levels - 1,
                        nrelaxc); /* nrelax on coarsest level */
    }
    core.SetCFactor(-1, cfactor);
    core.SetAccessLevel(access_level);
    core.SetPrintLevel(print_level);
    core.SetMaxIter(maxiter);
    core.SetAbsTol(tol);

    /* Parallel-in-time TriMGRIT simulation */
    start = clock();
    core.Drive();
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

    MPI_Finalize();

    return (0);
}
