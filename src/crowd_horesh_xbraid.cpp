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

#define CROWD_HORESH_LIBRARY
#include "crowd_horesh.cpp"
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

// #include "lapacke.h"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

class BraidVector
{
    public:
    Vector dm, drho, dlambda;

    BraidVector(Vector dm_, Vector drho_, Vector dlambda_) : dm(dm_), drho(drho_), dlambda(dlambda_) { }

    virtual ~BraidVector() { };
};

class MyBraidApp : public TriBraidApp {
protected:
    // BraidApp defines tstart, tstop, ntime and comm_t
public:
    int ntime, mspace, ilower, iupper, npoints, myid;
    double dx;
    FILE *ufile, *vfile, *wfile;
    std::vector<Vector> m, rho, lambda, q;
    Sparse K;

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

    const Vector apply_Phi(Vector &u, Vector &v,
                 const double t, const double dt);
};

MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_, double tstop_, int ntime_) : TriBraidApp(comm_t_, tstart_, tstop_, ntime_) {

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
    BraidVector *f = (BraidVector *) f_;
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

    Sparse X ((mspace + 1), mspace);
    
    X.insert(0, 0) = 0.25;
    for (int i = 1; i < mspace + 1; i++) {
        X.insert(i, i) = 0.25;
        X.insert(i - 1, i) = 0.25;
    }
    
    Sparse Pi = 2.0 * (X * (rho[index].cwiseInverse() + uright.)).asDiagonal();
    Sparse Qi = ...;

    return 0;
}

/*------------------------------------*/

/* Solve A(u) = f */

int MyBraidApp::TriSolve(braid_Vector uleft_, braid_Vector uright_,
                braid_Vector fleft_, braid_Vector fright_, braid_Vector f_,
                braid_Vector u_, BraidTriStatus &status) {
    BraidVector *uleft = (BraidVector *) uleft_;
    BraidVector *uright = (BraidVector *) uright_;
    BraidVector *f = (BraidVector *) f_;
    BraidVector *u = (BraidVector *) u_;
    BraidVector *fleft = (BraidVector *) fleft_;
    BraidVector *fright = (BraidVector *) fright_;

    double t, tprev, tnext, dt;

    /* Get the time-step size */
    status.GetTriT(&t, &tprev, &tnext);
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }

    // TODO

    /* no refinement */
    status.SetRFactor(1);

    return 0;
}

/*------------------------------------*/

/* This is only called from level 0 */

int MyBraidApp::Init(double t, braid_Vector *u_ptr_) {
    BraidVector **u_ptr = (BraidVector **) u_ptr_;

    // TODO

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Clone(braid_Vector u_, braid_Vector *v_ptr_) {
    BraidVector *u = (BraidVector *) u_;
    BraidVector **v_ptr = (BraidVector **) v_ptr_;

    // TODO

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

    // TODO

    return 0;
}

/*------------------------------------*/

int MyBraidApp::SpatialNorm(braid_Vector u_, double *norm_ptr) {
    BraidVector *u = (BraidVector *) u_;

    // TODO

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

        // TODO

        // fprintf(ufile, "%05d: ", index);
        // fprintf(vfile, "%05d: ", index);
        // fprintf(wfile, "%05d: ", index);
        for (j = 0; j < (mspace - 1); j++) {
            // fprintf(ufile, "% 1.14e, ", (*u->u)[j]);
            // fprintf(vfile, "% 1.14e, ", (*u->v)[j]);
            // fprintf(wfile, "% 1.14e, ", (*u->w)[j]);
        }
        // fprintf(ufile, "% 1.14e\n", (*u->u)[j]);
        // fprintf(vfile, "% 1.14e\n", (*u->v)[j]);
        // fprintf(wfile, "% 1.14e\n", (*u->w)[j]);
    }

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufSize(int *size_ptr, BraidBufferStatus &bstatus) {
    // *size_ptr = 3 * mspace * sizeof(double);
    // TODO
    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufPack(braid_Vector u_, void *buffer,
               BraidBufferStatus &bstatus) {
    BraidVector *u = (BraidVector *) u_;
    double *dbuffer = (double *) buffer;

    // TODO

    // bstatus.SetSize(3 * mspace * sizeof(double));

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufUnpack(void *buffer, braid_Vector *u_ptr_,
                 BraidBufferStatus &bstatus) {
    BraidVector **u_ptr = (BraidVector **) u_ptr_;

    // TODO
    bstatus.SetSize(3 * mspace * sizeof(double));

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
    // start = clock();
    core.Drive();
    // end = clock();
    mspace = 100;
    ntime = 150;
    time = 1.0;
    int iters = 15;

    double d_time = time / ntime;

    int START_LOOP_ITER = 1;
        
    app.m = std::vector<Vector>();
    for (int i = 0; i < ntime; i++) {
        Vector m_val (mspace + 1);
        m_val.setConstant(0.0);
        app.m.push_back(m_val);
    }

    app.rho = std::vector<Vector>();
    for (int i = 0; i < ntime + 1; i++) {
        Vector rho_val (mspace);
        rho_val.setConstant(0.5);
        app.rho.push_back(rho_val);
    }

    app.lambda = std::vector<Vector>();
    for (int i = 0; i < ntime + 2; i++) {
        Vector lambda_val (mspace);
        lambda_val.setConstant(0.1);
        app.lambda.push_back(lambda_val);
    }
    
    app.q = std::vector<Vector>(ntime + 2);

    double accumulator = 0.0;
    Vector q_val (mspace);
    for (int i = 0; i < mspace; i++) {
        q_val[i] = initial_condition(accumulator);
        accumulator += 1.0 / ((double)mspace - 1.0);
    }
    app.q[0] = q_val / d_time;

    accumulator = 0.0;
    q_val = Vector(mspace);
    for (int i = 0; i < mspace; i++) {
        q_val[i] = final_condition(accumulator);
        accumulator += 1.0 / ((double)mspace - 1.0);
    }
    app.q[app.q.size() - 1] = q_val / d_time;

    app.K = Vector(mspace, mspace + 1);

    for (int i = 0; i < mspace; i++) {
        app.K.insert(i, i) = -1.0 / dx;
        app.K.insert(i, i + 1) = 1.0 / dx;
    }

    core.Drive();

    /* Print runtime to file (for runtime comparisons)*/
    // time = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Total Run Time: %f s \n", time);

    MPI_Finalize();

    return 0;
}
