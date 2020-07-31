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
#define TAG_WORKER_RESULT 42
#define TYPE_RHO 2
#define TYPE_LAMBDA 3
#define TYPE_M 4

#define DM_LEN_SPACE (mspace+1)
#define DRHO_LEN_SPACE (mspace)
#define DLAMBDA_LEN_SPACE (mspace)
#define Q_LEN_SPACE (mspace)

#define DM_LEN_TIME (ntime)
#define DRHO_LEN_TIME (ntime + 1)
#define DLAMBDA_LEN_TIME (ntime + 2)
#define Q_LEN_TIME (ntime + 2)

// #include "lapacke.h"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

class BraidVector
{
    public:
    int index = -1;
    /// dm is length mspace + 2,
    /// drho is length mspace + 1,
    /// dlambda is length mspace
    Vector dm, drho, dlambda;
    /// Correct:
    /// dm is length mspace + 1,     (in time: ntime)
    /// drho is length mspace,       (in time: ntime + 1)
    /// dlambda is length mspace     (in time: ntime + 2)

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
    Vector GlambdaL, GmL, GrhoL;

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
    
    // Compute \nabla_m_i
    const Vector compute_GwMi(int index);
    // Compute \nabla_\lambda_i
    const Vector compute_GwLi(int index);
    // Compute \nabla_\rho_{i - 1/2}
    const Vector compute_GwRhoi(int index);
};

const Vector MyBraidApp::compute_GwMi(int index) {
    Eigen::Map<Vector> res(this->GmL.data() + (index * DM_LEN_SPACE), DM_LEN_SPACE);

    Vector res_(res);
    return res_;
}

const Vector MyBraidApp::compute_GwLi(int index) {
    Eigen::Map<Vector> res(this->GlambdaL.data() + (index * DLAMBDA_LEN_SPACE), DLAMBDA_LEN_SPACE);

    Vector res_(res);
    return res_;
}

const Vector MyBraidApp::compute_GwRhoi(int index) {
    Eigen::Map<Vector> res(this->GrhoL.data() + (index * DRHO_LEN_SPACE), DRHO_LEN_SPACE);

    Vector res_(res);
    return res_;
}

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
    int level;
    
    int index, final_index;

    status.GetTIndex(&index);
    status.GetNTPoints(&final_index);
    final_index -= 1;

    status.GetTriT(&t, &tprev, &tnext);
    status.GetLevel(&level);

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

    Vector mi = m[index];
    Vector mim1(mi.size());

    Vector rhoi = rho[index];
    Vector rhoip1(rhoi.size());

    if ((unsigned long) index == rho.size() - 1) {
        rhoip1.setZero();
    } else {
        rhoip1 = rho[index + 1];
    }
    
    // Compute Qi and Pi
    Vector tmp1 = rhoi.cwiseInverse() + rhoip1.cwiseInverse();
    auto Pi_ = 2.0 * ((X * tmp1).eval()).asDiagonal();
    Vector tmp2 = mi * mi + mim1 * mim1;
    Vector tmp3 = rhoi.cwiseProduct(rhoi.cwiseProduct(rhoi)).cwiseInverse();
    Vector tmp4 = X.transpose() * tmp2;
    Vector Qi_ (tmp3.size());
    Qi_.setConstant(2.0);
    Qi_ *= tmp3.asDiagonal();
    Qi_ *= tmp4.asDiagonal();
    auto Qi__ = Qi_.asDiagonal();
    Sparse Qi (Qi__);
    Sparse Pi (Pi_);

    Vector nabla_m_i = this->compute_GwMi(index);
    Vector nabla_rho_i = this->compute_GwMi(index);
    Vector nabla_lambda_i = this->compute_GwLi(index);

    r->dm = Pi * f->dm + K.transpose() * f->dlambda + nabla_m_i;
    r->drho = Qi * f->drho + uleft->dlambda / dt - f->dlambda / dt + nabla_rho_i;
    r->dlambda = K * f->dm + uright->drho / dt - f->drho / dt + nabla_lambda_i;

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
    
    int index, final_index;

    status.GetTIndex(&index);
    status.GetNTPoints(&final_index);
    final_index -= 1;

    Sparse X ((mspace + 1), mspace);
    
    X.insert(0, 0) = 0.25;
    for (int i = 1; i < mspace; i++) {
        X.insert(i, i) = 0.25;
        X.insert(i - 1, i) = 0.25;
    }

    Vector mi = m[index];
    Vector mim1(mi.size());

    assert(index != 0);

    if (index == 0) {
        mim1.setZero();
    } else {
        mim1 = m[index - 1];
    }

    Vector rhoi = rho[index];
    // assert(rhoi.size() == DRHO_LEN_SPACE);
    Vector rhoip1(rhoi.size());

    if (mi.size() != mim1.size()) {
        printf("What: %ld, %ld, index=%d\n", mi.size(), mim1.size(), index);
        assert(false);
    }

    if ((unsigned long) index == rho.size() - 1) {
        rhoip1.setZero();
    } else {
        rhoip1 = rho[index + 1];
    }
    
    // Compute Qi and Pi
    Vector tmp1 = rhoi.cwiseInverse() + rhoip1.cwiseInverse();
    Vector tmp1_ = X * tmp1;
    auto Pi_ = 2.0 * (tmp1_).asDiagonal();
    assert(mi.cwiseProduct(mi).size() == mim1.cwiseProduct(mim1).size());
    Vector tmp2 = mi.cwiseProduct(mi) + mim1.cwiseProduct(mim1);
    Vector tmp3 = rhoi.cwiseProduct(rhoi.cwiseProduct(rhoi)).cwiseInverse();
    Vector tmp4 = X.transpose() * tmp2;
    Vector Qi_ (tmp3.size());
    Qi_.setConstant(2.0);
    Qi_ = Qi_.cwiseProduct(tmp3).cwiseProduct(tmp4);
    auto Qi__ = Qi_.asDiagonal();
    Sparse Qi (Qi__);
    Sparse Pi (Pi_);

    if (index == 1) {
        // Compute delta_rho_1
        // Equation 12
        Vector nabla_l_0 = this->compute_GwLi(0);
        u->drho = dt * -nabla_l_0;
        
        // Compute delta_lambda_1
        // Equation 13
        Vector nabla_rho_1 = this->compute_GwRhoi(1);
        // uleft should never be null here!
        assert(uleft != nullptr);
        Vector delta_lambda_0 = uleft->dlambda;
        Vector delta_rho_1 = u->drho;
        // Solve for delta_lambda_1
        u->dlambda = dt * (-nabla_rho_1 - delta_lambda_0 / dt);
        u->dlambda = dt * Qi * delta_rho_1;
    
        // Compute delta_m
        // Equation 9
        Vector nabla_m_i = this->compute_GwMi(index);
        // Setup Ax = b system
        Vector b = -nabla_m_i - K.transpose() * u->dlambda;
        
        // Solve Pi x = b
        u->dm = Pi.diagonal().cwiseInverse().cwiseProduct(b);
    } else {
        // Compute delta_lambda_i
        // Equation 24
        Vector nabla_m_i = this->compute_GwMi(index);
        Vector nabla_lambda_i = this->compute_GwLi(index);
        Vector nabla_rho_i = this->compute_GwRhoi(index);
        Vector delta_lambda_im1;
        if (uleft == nullptr) {
            delta_lambda_im1 = Vector(u->dlambda.size());
            delta_lambda_im1.setZero();
        } else {
            delta_lambda_im1 = uleft->dlambda;
        }
        Vector delta_rho_ip1;
        if (uright == nullptr) {
            delta_rho_ip1 = Vector(u->drho.size());
            delta_rho_ip1.setZero();
        } else {
            delta_rho_ip1 = uright->drho;
        }
        Sparse I (DLAMBDA_LEN_SPACE, DLAMBDA_LEN_SPACE);
        I.setIdentity();
        Sparse A = -dt * Qi * K * Pi.cwiseInverse() * K.transpose();
        A -= I / dt;
        Vector b = dt * Qi * K * Pi.cwiseInverse() * nabla_m_i - Qi * delta_rho_ip1 - dt * Qi * nabla_lambda_i - delta_lambda_im1 / dt - nabla_rho_i;
        u->dlambda = A.cwiseInverse() * b;

        // Compute delta_m
        // Equation 9
        if (index != 0 && index != final_index) { // Otherwise delta_m doesn't matter
            Vector nabla_m_i = this->compute_GwMi(index);
            // Setup Ax = b system
            Vector b = -nabla_m_i - K.transpose() * u->dlambda;
            
            // Solve Pi x = b
            u->dm = Pi.cwiseInverse() * b;
        }

        // Compute delta_rho_i
        if (index != 0) {
            // Equation 11
            u->drho = Qi.cwiseInverse() * (-nabla_rho_i - delta_lambda_im1 / dt + u->dlambda / dt);
        }
    }

    /* no refinement */
    status.SetRFactor(1);

    return 0;
}

/*------------------------------------*/

/* This is only called from level 0 */

int MyBraidApp::Init(double t, braid_Vector *u_ptr_) {
    BraidVector **u_ptr = (BraidVector **) u_ptr_;

    Vector dm(DM_LEN_SPACE);
    Vector drho(DRHO_LEN_SPACE);
    Vector dlambda(DLAMBDA_LEN_SPACE);

    *u_ptr = new BraidVector(dm, drho, dlambda);

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Clone(braid_Vector u_, braid_Vector *v_ptr_) {
    BraidVector *u = (BraidVector *) u_;
    BraidVector **v_ptr = (BraidVector **) v_ptr_;
    
    // TODO: maybe?:
    // delete *v_ptr;

    *v_ptr = new BraidVector(u->dm, u->drho, u->dlambda);
    (*v_ptr)->index = u->index;

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

    y->dm = alpha * x->dm + beta * y->dm;
    y->drho = alpha * x->drho + beta * y->drho;
    y->dlambda = alpha * x->dlambda + beta * y->dlambda;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::SpatialNorm(braid_Vector u_, double *norm_ptr) {
    BraidVector *u = (BraidVector *) u_;

    *norm_ptr = sqrt(u->dlambda.squaredNorm() + u->dm.squaredNorm() + u->drho.squaredNorm());

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Access(braid_Vector u_, BraidAccessStatus &astatus) {
    BraidVector *u = (BraidVector *) u_;

    int done, index;

    astatus.GetDone(&done);
    if (done) {
        astatus.GetTIndex(&index);

        MPI_Request send_request;
        int message_size = sizeof(int) * 2 + sizeof(double) * (mspace + 2);
        
        // type is TYPE_RHO
        {
        char *buffer = (char *) calloc(message_size, 1);
        int *ibuffer = (int *) buffer;
        *(ibuffer++) = TYPE_RHO;
        *(ibuffer++) = index;
        double *dbuffer = (double *) ibuffer;
        for (int i = 0; i < DRHO_LEN_SPACE; i++) {
            dbuffer[i] = u->drho[i];
        }
        MPI_Isend((void *) buffer, message_size, MPI_BYTE, 0, TAG_WORKER_RESULT, MPI_COMM_WORLD, &send_request);
        free(buffer);
        }
        // type is TYPE_LAMBDA
        {
        char *buffer = (char *) calloc(message_size, 1);
        int *ibuffer = (int *) buffer;
        *(ibuffer++) = TYPE_LAMBDA;
        *(ibuffer++) = index;
        double *dbuffer = (double *) ibuffer;
        for (int i = 0; i < DLAMBDA_LEN_SPACE; i++) {
            dbuffer[i] = u->dlambda[i];
        }
        MPI_Isend((void *) buffer, message_size, MPI_BYTE, 0, TAG_WORKER_RESULT, MPI_COMM_WORLD, &send_request);
        free(buffer);
        }
        // type is TYPE_M
        {
        char *buffer = (char *) calloc(message_size, 1);
        int *ibuffer = (int *) buffer;
        *(ibuffer++) = TYPE_M;
        *(ibuffer++) = index;
        double *dbuffer = (double *) ibuffer;
        for (int i = 0; i < DM_LEN_SPACE; i++) {
            dbuffer[i] = u->dm[i];
        }
        MPI_Isend((void *) buffer, message_size, MPI_BYTE, 0, TAG_WORKER_RESULT, MPI_COMM_WORLD, &send_request);
        free(buffer);
        }
    }

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufSize(int *size_ptr, BraidBufferStatus &bstatus) {
    // sizeof(index) + sizeof(drho) + sizeof(dlambda) + sizeof(dm)
    *size_ptr = sizeof(int) + DRHO_LEN_SPACE * sizeof(double) + DLAMBDA_LEN_SPACE * sizeof(double) + DM_LEN_SPACE * sizeof(double);
    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufPack(braid_Vector u_, void *buffer_,
               BraidBufferStatus &bstatus) {
    BraidVector *u = (BraidVector *) u_;
    int *buffer = (int *) buffer_;
    
    *(buffer++) = u->index;

    double *dbuffer = (double *) buffer;

    for (int i = 0; i < DRHO_LEN_SPACE; i++, dbuffer++) {
        *dbuffer = u->drho[i];
    }

    for (int i = 0; i < DLAMBDA_LEN_SPACE; i++, dbuffer++) {
        *dbuffer = u->dlambda[i];
    }

    for (int i = 0; i < DM_LEN_SPACE; i++, dbuffer++) {
        *dbuffer = u->dm[i];
    }

    bstatus.SetSize(sizeof(int) + DRHO_LEN_SPACE * sizeof(double) + DLAMBDA_LEN_SPACE * sizeof(double) + DM_LEN_SPACE * sizeof(double));

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufUnpack(void *buffer_, braid_Vector *u_ptr_,
                 BraidBufferStatus &bstatus) {
    BraidVector **u_ptr = (BraidVector **) u_ptr_;
    
    int *buffer = (int *) buffer_;

    (*u_ptr)->index = *(buffer++);

    double *dbuffer = (double *) buffer;

    for (int i = 0; i < DRHO_LEN_SPACE; i++, dbuffer++) {
        (*u_ptr)->drho[i] = *dbuffer;
    }

    for (int i = 0; i < DLAMBDA_LEN_SPACE; i++, dbuffer++) {
        (*u_ptr)->dlambda[i] = *dbuffer;
    }

    for (int i = 0; i < DM_LEN_SPACE; i++, dbuffer++) {
        (*u_ptr)->dm[i] = *dbuffer;
    }

    bstatus.SetSize(sizeof(int) + DRHO_LEN_SPACE * sizeof(double) + DLAMBDA_LEN_SPACE * sizeof(double) + DM_LEN_SPACE * sizeof(double));

    return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    double tstart, tstop, dx;
    int rank, ntime, mspace;
    int max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
    int access_level, print_level;
    double tol;
    double time;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Define space domain. Space domain is between 0 and 1, mspace defines the
     * number of steps */
    mspace = 100;
    ntime = 150;

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
    time = 1.0;
    int iters = 15;

    double d_time = time / ntime;


    app.m = std::vector<Vector>();
    for (int i = 0; i < DM_LEN_TIME; i++) {
        Vector m_val (DM_LEN_SPACE);
        m_val.setConstant(0.0);
        app.m.push_back(m_val);
    }

    app.rho = std::vector<Vector>();
    for (int i = 0; i < DRHO_LEN_TIME; i++) {
        Vector rho_val (DRHO_LEN_SPACE);
        rho_val.setConstant(0.5);
        app.rho.push_back(rho_val);
    }

    app.lambda = std::vector<Vector>();
    for (int i = 0; i < DLAMBDA_LEN_TIME; i++) {
        Vector lambda_val (DLAMBDA_LEN_SPACE);
        lambda_val.setConstant(0.1);
        app.lambda.push_back(lambda_val);
    }
    
    app.q = std::vector<Vector>(Q_LEN_TIME);

    double accumulator = 0.0;
    Vector q_val (Q_LEN_SPACE);
    for (int i = 0; i < mspace; i++) {
        q_val[i] = initial_condition(accumulator);
        accumulator += 1.0 / ((double)mspace - 1.0);
    }
    app.q[0] = q_val / d_time;

    for (int i = 1; i < ntime + 1; i++) {
        app.q[i] = q_val * 0.0;
    }

    accumulator = 0.0;
    q_val = Vector(mspace);
    for (int i = 0; i < mspace; i++) {
        q_val[i] = final_condition(accumulator);
        accumulator += 1.0 / ((double)mspace - 1.0);
    }
    app.q[app.q.size() - 1] = q_val / d_time;

    app.K = Sparse(mspace, mspace + 1);

    Sparse D1 = get_derivative_matrix_space(mspace, ntime, d_time);
    Sparse D2 = get_derivative_matrix_time(mspace, ntime, d_time);
    Sparse D = joinlr(D1, D2);
    Sparse As = get_As(mspace, ntime);
    Sparse At = get_At(mspace, ntime);

    for (int i = 0; i < mspace; i++) {
        app.K.insert(i, i) = -1.0 / dx;
        app.K.insert(i, i + 1) = 1.0 / dx;
    }

    Vector q_long (app.q.size() * app.q[0].size());
    for (unsigned long i = 0; i < app.q.size(); i++) {
        for (int j = 0; j < app.q[0].size(); j++) {
            double val = app.q[i][j];
            q_long[i * app.q[0].size() + j] = val;
        }
    }

    for (int i = 0; i < iters; i++) {
        Vector rho_long (app.rho.size() * app.rho[0].size());
        for (unsigned long i = 0; i < app.rho.size(); i++) {
            for (int j = 0; j < app.rho[0].size(); j++) {
                rho_long[i * app.rho[0].size() + j] = app.rho[i][j];
            }
        }
        Vector m_long (app.m.size() * app.m[0].size());
        for (unsigned long i = 0; i < app.m.size(); i++) {
            for (int j = 0; j < app.m[0].size(); j++) {
                m_long[i * app.m[0].size() + j] = app.m[i][j];
            }
        }
        Vector lambda_long (app.lambda.size() * app.lambda[0].size());
        for (unsigned long i = 0; i < app.lambda.size(); i++) {
            for (int j = 0; j < app.lambda[0].size(); j++) {
                lambda_long[i * app.lambda[0].size() + j] = app.lambda[i][j];
            }
        }
        
        app.GlambdaL = get_GlambdaL(m_long, rho_long, q_long, D);
        Vector GmL1 = 2.0 * m_long.asDiagonal() * As.transpose() * At * rho_long.cwiseInverse();
        Vector GmL2 = D1.transpose() * lambda_long;
    
        app.GmL = GmL1 + GmL2;
        app.GrhoL =
            rho_long.cwiseProduct(rho_long).eval().cwiseInverse().eval().asDiagonal() *
                (-At.transpose()) * As * m_long.cwiseProduct(m_long) +
            D2.transpose() * lambda_long;
        

        core.Drive();

        if (rank == 0) {
            // Main processor
            // Receive from access the results of workers' computation
            std::vector<Vector> drho(DRHO_LEN_TIME);
            std::vector<Vector> dlambda(DLAMBDA_LEN_TIME);
            std::vector<Vector> dm(DM_LEN_TIME);

            int num_received = 0;

            while (num_received < DRHO_LEN_TIME * DLAMBDA_LEN_TIME * DM_LEN_TIME) {
                int type, index;
                int message_size = sizeof(int) * 2 + sizeof(double) * (mspace + 2);
                char *buffer = (char *) malloc(message_size);
                MPI_Recv((void *) buffer, message_size, MPI_BYTE, MPI_ANY_SOURCE, TAG_WORKER_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int *buffer_ = (int *) buffer;
                type = *(buffer++);
                index = *(buffer++);
                double *dbuffer = (double *) buffer_;
                if (type == TYPE_RHO) {
                    Vector drho_val(DRHO_LEN_SPACE);
                    for (int i = 0; i < DRHO_LEN_SPACE; i++) {
                        drho_val[i] = dbuffer[i];
                    }
                    drho[index] = drho_val;
                } else if (type == TYPE_LAMBDA) {
                    Vector dlambda_val(DLAMBDA_LEN_SPACE);
                    for (int i = 0; i < DLAMBDA_LEN_SPACE; i++) {
                        dlambda_val[i] = dbuffer[i];
                    }
                    dlambda[index] = dlambda_val;
                } else if (type == TYPE_M) {
                    Vector dm_val(DM_LEN_SPACE);
                    for (int i = 0; i < DM_LEN_SPACE; i++) {
                        dm_val[i] = dbuffer[i];
                    }
                    dm[index] = dm_val;
                }
                num_received++;
                free(buffer);
            }
            // Compute global new m, rho, and lambda using line search etc.
            Vector drho_long (drho.size() * drho[0].size());
            for (unsigned long i = 0; i < drho.size(); i++) {
                for (int j = 0; j < drho[0].size(); j++) {
                    drho_long[i * drho[0].size() + j] = drho[i][j];
                }
            }
            Vector dm_long (dm.size() * dm[0].size());
            for (unsigned long i = 0; i < dm.size(); i++) {
                for (int j = 0; j < dm[0].size(); j++) {
                    dm_long[i * dm[0].size() + j] = dm[i][j];
                }
            }
            Vector dlambda_long (dlambda.size() * dlambda[0].size());
            for (unsigned long i = 0; i < dlambda.size(); i++) {
                for (int j = 0; j < dlambda[0].size(); j++) {
                    dlambda_long[i * dlambda[0].size() + j] = dlambda[i][j];
                }
            }
            double alpha = line_search(dm_long, drho_long, dlambda_long, m_long, rho_long, lambda_long, As, At, D1, D2, q_long, D);
            
            for (unsigned long i = 0; i < app.m.size(); i++) {
                app.m[i] += alpha * dm[i];
            }
            for (unsigned long i = 0; i < app.m.size(); i++) {
                app.rho[i] += alpha * drho[i];
            }
            for (unsigned long i = 0; i < app.m.size(); i++) {
                app.lambda[i] += alpha * dlambda[i];
            }
            // Send global picture of m, rho, and lambda
            int message_length = sizeof(double) * (app.m.size() * app.m[0].size() + app.rho.size() * app.rho[0].size() + app.lambda.size() * app.lambda[0].size());
            double *buffer = (double *) malloc(message_length * sizeof(double));
            double *buffer_ptr = buffer;
            for (Vector &m : app.m) {
                for (double &val : m) {
                    *(buffer_ptr++) = val;
                }
            }
            for (Vector &rho : app.rho) {
                for (double &val : rho) {
                    *(buffer_ptr++) = val;
                }
            }
            for (Vector &lambda : app.lambda) {
                for (double &val : lambda) {
                    *(buffer_ptr++) = val;
                }
            }
            MPI_Bcast(buffer, message_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            free(buffer);
        } else {
            // Worker
            // access should have already sent their results to the main processor
            // Receive a global picture of m, rho, and lambda
            int message_length = sizeof(double) * (app.m.size() * app.m[0].size() + app.rho.size() * app.rho[0].size() + app.lambda.size() * app.lambda[0].size());
            double *buffer = (double *) malloc(message_length * sizeof(double));
            MPI_Bcast(buffer, message_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            double *buffer_ptr = buffer;
            for (Vector &m : app.m) {
                for (double &val : m) {
                    val = *(buffer_ptr++);
                }
            }
            for (Vector &rho : app.rho) {
                for (double &val : rho) {
                    val = *(buffer_ptr++);
                }
            }
            for (Vector &lambda : app.lambda) {
                for (double &val : lambda) {
                    val = *(buffer_ptr++);
                }
            }
            free(buffer);
        }
    }

    // TODO: 
    // rank == 0 should print some results

    /* Print runtime to file (for runtime comparisons)*/
    // time = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Total Run Time: %f s \n", time);
    
    printf("Hit the end. TODO: Print out actual results\n");

    MPI_Finalize();

    return 0;
}
