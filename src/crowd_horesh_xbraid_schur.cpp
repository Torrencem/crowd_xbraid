/*HEADER**********************************************************************
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

// Crowd model problem
//
// This program simulates an optimization problem of modelling the movement of a
// crowd from one configuration to another. The math behind the code used here
// can be found in (TODO: Reference to haber and horesh)

#include "horesh_utils.cpp"
#include "initial_final.cpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "braid_test.h"
#include "tribraid.hpp"
#define TAG_WORKER_RESULT 42

#define DM_LEN_SPACE (mspace + 1)
#define DRHO_LEN_SPACE (mspace)
#define DLAMBDA_LEN_SPACE (mspace)
#define Q_LEN_SPACE (mspace)

#define DM_LEN_TIME (ntime)
#define DRHO_LEN_TIME (ntime + 1)
#define DLAMBDA_LEN_TIME (ntime + 2)
#define Q_LEN_TIME (ntime + 2)

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

class BraidVector {
  public:
    Vector dlambda;

    BraidVector(Vector &dlambda_) : dlambda(dlambda_) {}

    BraidVector(int mspace) {
        this->dlambda = Vector(DLAMBDA_LEN_SPACE);
        this->dlambda.setZero();
    }

    virtual ~BraidVector(){};
};

class MyBraidApp : public TriBraidApp {
  protected:
    // BraidApp defines tstart, tstop, ntime and comm_t
  public:
    int ntime, mspace, ilower, iupper, npoints, myid, cfactor;
    double dx;
    double normcoeff;
    std::vector<Vector> m, rho, lambda, q;
    Sparse K, X;
    Vector GlambdaL, GmL, GrhoL, RHS;

    // Main constructor
    MyBraidApp(MPI_Comm comm_t__, int rank_, double tstart_ = 0.0,
               double tstop_ = 1.0, int ntime_ = 100);

    virtual ~MyBraidApp(){};

    // Define the Vraid Wrapper routines
    // Note: braid_Vector is the type BraidVector*
    virtual int Clone(braid_Vector u_, braid_Vector *v_ptr) override;

    virtual int Init(double t, braid_Vector *v_ptr) override;

    virtual int Free(braid_Vector u_) override;

    virtual int Sum(double alpha, braid_Vector x_, double beta,
                    braid_Vector y_) override;

    virtual int SpatialNorm(braid_Vector u_, double *norm_ptr) override;

    virtual int Access(braid_Vector u_, BraidAccessStatus &astatus) override;

    virtual int TriResidual(braid_Vector uleft_, braid_Vector uright_,
                            braid_Vector f_, braid_Vector r_,
                            BraidTriStatus &pstatus) override;

    virtual int TriSolve(braid_Vector uleft_, braid_Vector uright_,
                         braid_Vector fleft_, braid_Vector fright_,
                         braid_Vector f_, braid_Vector u_,
                         BraidTriStatus &pstatus) override;

    virtual int BufSize(braid_Int *size_ptr,
                        BraidBufferStatus &bstatus) override;

    virtual int BufPack(braid_Vector u_, void *buffer,
                        BraidBufferStatus &bstatus) override;

    virtual int BufUnpack(void *buffer, braid_Vector *u_ptr,
                          BraidBufferStatus &bstatus) override;

    const Sparse computeP(const int index);
    const Sparse computeQ(const int index);
    const Vector get_RHS(const int index);
    // Not needed
    virtual int Residual(braid_Vector u_, braid_Vector r_,
                         BraidStepStatus &pstatus) override {
        return 0;
    }

    const Vector apply_Phi(Vector &u, Vector &v, const double t,
                           const double dt);

    // Compute \nabla_m_i
    const Vector compute_GwMi(const int index);
    // Compute \nabla_\lambda_i
    const Vector compute_GwLi(const int index);
    // Compute \nabla_\rho_{i - 1/2}
    const Vector compute_GwRhoi(const int index);
};

const Vector MyBraidApp::get_RHS(int index) {
    return RHS(Eigen::seq(index * DLAMBDA_LEN_SPACE,
                          (index + 1) * DLAMBDA_LEN_SPACE - 1));
}

const Vector MyBraidApp::compute_GwMi(int index) {
    Eigen::Map<Vector> res(this->GmL.data() + (index * DM_LEN_SPACE),
                           DM_LEN_SPACE);

    Vector res_(res);
    return res_;
}

const Vector MyBraidApp::compute_GwLi(const int index) {
    Eigen::Map<Vector> res(this->GlambdaL.data() + (index * DLAMBDA_LEN_SPACE),
                           DLAMBDA_LEN_SPACE);

    Vector res_(res);
    return res_;
}

const Vector MyBraidApp::compute_GwRhoi(const int index) {
    Eigen::Map<Vector> res(this->GrhoL.data() + (index * DRHO_LEN_SPACE),
                           DRHO_LEN_SPACE);

    Vector res_(res);
    return res_;
}

MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_,
                       double tstop_, int ntime_)
    : TriBraidApp(comm_t_, tstart_, tstop_, ntime_) {}

//------------------------------------

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

// Compute the matrix P_i, which forms one chunk of the top left quadrant of
// \hat{A}.
const Sparse MyBraidApp::computeP(int index) {
    assert(index < DRHO_LEN_TIME);
    Vector tmp1 = rho[index].cwiseInverse() + rho[index + 1].cwiseInverse();
    Vector tmp1_ = X * tmp1;
    auto Pi_ = 2.0 * (tmp1_).asDiagonal();
    return Sparse(Pi_);
}

// Compute the matrix Q_i, which forms one chunk of the bottom right quadrant of
// \hat{A}.
const Sparse MyBraidApp::computeQ(int index) {

    Vector mi(DM_LEN_SPACE);

    if (index == DM_LEN_TIME) {
        mi.setZero();
    } else {
        mi = m[index];
    }

    Vector mim1(DM_LEN_SPACE);

    if (index == 0) {
        mim1.setZero();
    } else {
        mim1 = m[index - 1];
    }

    Vector tmp2 = mi.cwiseProduct(mi) + mim1.cwiseProduct(mim1);
    Vector rhoi = rho[index];
    Vector tmp3 = rhoi.cwiseProduct(rhoi.cwiseProduct(rhoi)).cwiseInverse();
    Vector tmp4 = X.transpose() * tmp2;
    Vector Qi_(tmp3.size());
    Qi_.setConstant(2.0);
    Qi_ = Qi_.cwiseProduct(tmp3).cwiseProduct(tmp4);
    auto Qi__ = Qi_.asDiagonal();
    Sparse normalizer(Qi_.size(), Qi_.size());
    normalizer.setIdentity();
    return Sparse(Qi__) + normcoeff * normalizer;
}

// Compute r = A(r) - f

int MyBraidApp::TriResidual(braid_Vector uleft_, braid_Vector uright_,
                            braid_Vector f_, braid_Vector r_,
                            BraidTriStatus &status) {
    BraidVector *uleft = (BraidVector *)uleft_;
    BraidVector *uright = (BraidVector *)uright_;
    BraidVector *f = (BraidVector *)f_;
    BraidVector *r = (BraidVector *)r_;

    double t, tprev, tnext, dt;

    int index, final_index;

    // Get the time-step size and the ``true" index of the vector in question,
    // compensating for coarsening.
    status.GetTIndex(&index);
    status.GetNTPoints(&final_index);
    status.GetTriT(&t, &tprev, &tnext);
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }
    int level;
    status.GetLevel(&level);
    int factor = pow(cfactor, level);
    index = index * factor;
    dt = factor / (factor / dt -
                   1); // Since we are calculating dlambda, we have one more
                       // TriMGRIT time points than we do real time points.

    if (uleft == nullptr) {
        Sparse Q0 = invertDiagonal(computeQ(0)) / (dt * dt);
        r->dlambda = Q0 * r->dlambda - Q0 * uright->dlambda;
    } else if (uright == nullptr) {
        Sparse Qn = invertDiagonal(computeQ(index - factor)) / (dt * dt);
        r->dlambda = -Qn * uleft->dlambda + Qn * r->dlambda;
    } else {
        Sparse Qi = invertDiagonal(computeQ(index)) / (dt * dt);
        Sparse Qim1 = invertDiagonal(computeQ(index - factor)) / (dt * dt);
        Sparse KPiKt =
            K * invertDiagonal(computeP(index - factor)) * K.transpose();
        r->dlambda = -Qim1 * uleft->dlambda + (Qi + Qim1 + KPiKt) * r->dlambda -
                     Qi * uright->dlambda;
    }

    Vector RHS = get_RHS(index);

    r->dlambda = r->dlambda - RHS;
    if (f != nullptr)
        r->dlambda = r->dlambda - f->dlambda;

    return 0;
}

//------------------------------------

// Solve A(u) = f

int MyBraidApp::TriSolve(braid_Vector uleft_, braid_Vector uright_,
                         braid_Vector fleft_, braid_Vector fright_,
                         braid_Vector f_, braid_Vector u_,
                         BraidTriStatus &status) {
    BraidVector *u = (BraidVector *)u_;
    BraidVector r = *u;

    // Calculate the residual in the equation and store it in r.
    TriResidual(uleft_, uright_, f_, (braid_Vector)&r, status);

    // See TriResidual for explanation of this block.
    double t, tprev, tnext, dt;

    status.GetTriT(&t, &tprev, &tnext);
    if (t < tnext) {
        dt = tnext - t;
    } else {
        dt = t - tprev;
    }

    int index, final_index;
    status.GetTIndex(&index);
    status.GetNTPoints(&final_index);
    int level;
    status.GetLevel(&level);
    int factor = pow(cfactor, level);
    index = index * factor;
    dt = factor / (factor / dt - 1);

    Sparse C;

    if (uleft_ == nullptr) {
        C = invertDiagonal(computeQ(0)) / (dt * dt);
    } else if (uright_ == nullptr) {
        C = invertDiagonal(computeQ(index - factor)) / (dt * dt);
    } else {
        Sparse Qi = invertDiagonal(computeQ(index)) / (dt * dt);
        Sparse Qim1 = invertDiagonal(computeQ(index - factor)) / (dt * dt);
        Sparse KPiKt =
            K * invertDiagonal(computeP(index - factor)) * K.transpose();
        C = Qi + Qim1 + KPiKt;
    }

    // Since C is tridiagonal, we could use a more efficient solver.
    Eigen::BiCGSTAB<Sparse, Eigen::IncompleteLUT<double>> solver;
    solver.compute(C);
    u->dlambda = u->dlambda - solver.solve(r.dlambda);

    status.SetRFactor(1); // No refinement

    return 0;
}

//------------------------------------

// This is only called from level 0

int MyBraidApp::Init(double t, braid_Vector *u_ptr_) {
    BraidVector **u_ptr = (BraidVector **)u_ptr_;

    *u_ptr = new BraidVector(mspace);

    (*u_ptr)->dlambda.setConstant(0.1);

    return 0;
}

//------------------------------------

int MyBraidApp::Clone(braid_Vector u_, braid_Vector *v_ptr_) {
    BraidVector *u = (BraidVector *)u_;
    BraidVector **v_ptr = (BraidVector **)v_ptr_;

    // TODO: maybe?:
    // delete *v_ptr;

    *v_ptr = new BraidVector(u->dlambda);

    return 0;
}

//------------------------------------

int MyBraidApp::Free(braid_Vector u_) {
    BraidVector *u = (BraidVector *)u_;
    delete u;

    return 0;
}

//------------------------------------

int MyBraidApp::Sum(double alpha, braid_Vector x_, double beta,
                    braid_Vector y_) {
    BraidVector *y = (BraidVector *)y_;
    BraidVector *x = (BraidVector *)x_;

    y->dlambda = alpha * x->dlambda + beta * y->dlambda;

    return 0;
}

//------------------------------------

int MyBraidApp::SpatialNorm(braid_Vector u_, double *norm_ptr) {
    BraidVector *u = (BraidVector *)u_;

    *norm_ptr = sqrt(u->dlambda.squaredNorm());

    return 0;
}

//------------------------------------

int MyBraidApp::Access(braid_Vector u_, BraidAccessStatus &astatus) {
    BraidVector *u = (BraidVector *)u_;

    int done, index;

    astatus.GetDone(&done);
    if (done) {
        astatus.GetTIndex(&index);
        //        std::cout << "At index " << index << std::endl;
        //        std::cout << "dlambda is " << u->dlambda << std::endl;

        MPI_Request send_request;
        int message_size = sizeof(int) + sizeof(double) * DLAMBDA_LEN_SPACE;
        char *buffer = (char *)calloc(message_size, 1);
        int *ibuffer = (int *)buffer;
        *(ibuffer++) = index;
        double *dbuffer = (double *)ibuffer;
        for (int i = 0; i < DLAMBDA_LEN_SPACE; i++) {
            dbuffer[i] = u->dlambda[i];
        }
        MPI_Isend((void *)buffer, message_size, MPI_BYTE, 0, TAG_WORKER_RESULT,
                  MPI_COMM_WORLD, &send_request);
        free(buffer);
    }
    return 0;
}

//------------------------------------

int MyBraidApp::BufSize(int *size_ptr, BraidBufferStatus &bstatus) {
    *size_ptr = DLAMBDA_LEN_SPACE * sizeof(double);

    return 0;
}

//------------------------------------

int MyBraidApp::BufPack(braid_Vector u_, void *buffer_,
                        BraidBufferStatus &bstatus) {
    BraidVector *u = (BraidVector *)u_;

    double *dbuffer = (double *)buffer_;

    for (int i = 0; i < DLAMBDA_LEN_SPACE; i++) {
        dbuffer[i] = u->dlambda[i];
    }

    bstatus.SetSize(DLAMBDA_LEN_SPACE * sizeof(double));

    return 0;
}

//------------------------------------

int MyBraidApp::BufUnpack(void *buffer_, braid_Vector *u_ptr_,
                          BraidBufferStatus &bstatus) {
    BraidVector **u_ptr = (BraidVector **)u_ptr_;

    *u_ptr = new BraidVector(mspace);

    double *dbuffer = (double *)buffer_;

    for (int i = 0; i < DLAMBDA_LEN_SPACE; i++) {
        double val = dbuffer[i];
        (*u_ptr)->dlambda[i] = val;
    }

    bstatus.SetSize(DLAMBDA_LEN_SPACE * sizeof(double));

    return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    double tstart, tstop, dx;
    int rank, ntime, mspace;
    int max_levels, min_coarse, nrelax, nrelaxc, maxiter;
    int access_level, print_level;
    double tol;
    double time;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Define space domain. Space domain is between 0 and 1, mspace defines the
    // number of steps.
    mspace = 8;
    ntime = 8;

    // Define some Braid parameters
    max_levels = 2;
    min_coarse = 1;
    nrelax = 25;
    nrelaxc = 25;
    maxiter = 6000;
    tol = 1.0e-10;
    access_level = 1;
    print_level = 2;

    // Define the space step
    dx = (double)1 / (mspace + 1);

    // Define time domain and step
    tstart = 0.0; // Beginning of time domain
    tstop = 1.0;  // End of time domain

    // Set up the app structure
    // Since we're calculating dlambda, we use ntime + 1 time points.
    auto app = MyBraidApp(MPI_COMM_WORLD, rank, tstart, tstop, ntime + 1);
    app.myid = rank;
    app.ntime = ntime;
    app.mspace = mspace;
    app.dx = dx;
    app.cfactor = 2;
    app.normcoeff = 0.0001; // The normalization coefficient used to ensure that
                            // Q_i is invertible.  This will approach zero as
                            // the algorithm iterates.

    // Initialize XBraid

    auto core = BraidTriCore(MPI_COMM_WORLD, &app);

    // Set some XBraid(_Adjoint) parameters
    core.SetMaxLevels(max_levels);
    core.SetMinCoarse(min_coarse);
    core.SetNRelax(-1, nrelax);
    if (max_levels > 1) {
        core.SetNRelax(max_levels - 1, nrelaxc); // nrelax on coarsest level
    }
    core.SetCFactor(-1, app.cfactor);
    core.SetAccessLevel(access_level);
    core.SetPrintLevel(print_level);
    core.SetMaxIter(maxiter);
    core.SetAbsTol(tol);

    time = 1.0;
    int iters = 2;

    double d_time = time / ntime;

    // Set up initial values for m, rho, and lambda.  (TriMGRIT won't see these:
    // it calculates dm, drho, and dlambda.)
    app.m = std::vector<Vector>();
    for (int i = 0; i < DM_LEN_TIME; i++) {
        Vector m_val(DM_LEN_SPACE);
        m_val.setConstant(0.0);
        app.m.push_back(m_val);
    }

    app.rho = std::vector<Vector>();
    for (int i = 0; i < DRHO_LEN_TIME; i++) {
        Vector rho_val(DRHO_LEN_SPACE);
        rho_val.setConstant(0.6);
        app.rho.push_back(rho_val);
    }

    app.lambda = std::vector<Vector>();
    for (int i = 0; i < DLAMBDA_LEN_TIME; i++) {
        Vector lambda_val(DLAMBDA_LEN_SPACE);
        lambda_val.setConstant(0.1);
        app.lambda.push_back(lambda_val);
    }

    // Set up initial and final values for rho.
    app.q = std::vector<Vector>(Q_LEN_TIME);

    double accumulator = 0.0;
    Vector q_val(Q_LEN_SPACE);
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
        q_val[i] = -1.0 * final_condition(accumulator);
        accumulator += 1.0 / ((double)mspace - 1.0);
    }
    app.q[app.q.size() - 1] = q_val / d_time;

    Vector q_long(app.q.size() * app.q[0].size());
    for (unsigned long i = 0; i < app.q.size(); i++) {
        for (int j = 0; j < app.q[0].size(); j++) {
            double val = app.q[i][j];
            q_long[i * app.q[0].size() + j] = val;
        }
    }

    // Various constant matrices used in the discretization.
    app.K = Sparse(mspace, mspace + 1);
    app.X = Sparse(mspace + 1, mspace);

    Sparse D1 = get_derivative_matrix_space(mspace, ntime, d_time);
    Sparse D2 = get_derivative_matrix_time(mspace, ntime, d_time);
    Sparse D = joinlr(D1, D2);
    Sparse As = get_As(mspace, ntime); // Averaging matrix in space
    Sparse At = get_At(mspace, ntime); // Averaging matrix in time

    for (int i = 0; i < mspace; i++) {
        app.K.insert(i, i) = -1.0 / dx + 1;
        app.K.insert(i, i + 1) = 1.0 / dx - 1;
    }

    for (int i = 0; i < mspace; i++) {
        app.X.insert(i, i) = 0.25;
        app.X.insert(i + 1, i) = 0.25;
    }

    // Iterative TriMGRIT SQP
    for (int i_ = 0; i_ < iters; i_++) {

        // Full Eigen vectors for rho, m, and lambda.
        Vector rho_long(DRHO_LEN_SPACE * DRHO_LEN_TIME);
        for (int i = 0; i < DRHO_LEN_TIME; i++) {
            for (int j = 0; j < DRHO_LEN_SPACE; j++) {
                rho_long[i * DRHO_LEN_SPACE + j] = app.rho[i][j];
            }
        }
        Vector m_long(DM_LEN_SPACE * DM_LEN_TIME);
        for (int i = 0; i < DM_LEN_TIME; i++) {
            for (int j = 0; j < DM_LEN_SPACE; j++) {
                m_long[i * DM_LEN_SPACE + j] = app.m[i][j];
            }
        }
        Vector lambda_long(DLAMBDA_LEN_SPACE * DLAMBDA_LEN_TIME);
        for (int i = 0; i < DLAMBDA_LEN_TIME; i++) {
            for (int j = 0; j < DLAMBDA_LEN_SPACE; j++) {
                lambda_long[i * DLAMBDA_LEN_SPACE + j] = app.lambda[i][j];
            }
        }

        // Compute various gradients, culminating in the right-hand side of the
        // equation we wish to solve.
        app.GlambdaL = get_GlambdaL(m_long, rho_long, q_long, D);
        Vector GmL1 = 2.0 * m_long.asDiagonal() * As.transpose() * At *
                      rho_long.cwiseInverse();
        Vector GmL2 = D1.transpose() * lambda_long;

        app.GmL = GmL1 + GmL2;
        app.GrhoL = rho_long.cwiseProduct(rho_long)
                            .eval()
                            .cwiseInverse()
                            .eval()
                            .asDiagonal() *
                        (-At.transpose()) * As * m_long.cwiseProduct(m_long) +
                    D2.transpose() * lambda_long + app.normcoeff * rho_long;

        Vector RHS_m(m_long.size());
        for (int i = 0; i < DM_LEN_TIME; i++) {
            auto sequence =
                Eigen::seq(i * DM_LEN_SPACE, (i + 1) * DM_LEN_SPACE - 1);
            RHS_m(sequence)
                << invertDiagonal(app.computeP(i)) * app.GmL(sequence);
        }
        Vector RHS_rho(rho_long.size());
        for (int i = 0; i < DRHO_LEN_TIME; i++) {
            auto sequence =
                Eigen::seq(i * DRHO_LEN_SPACE, (i + 1) * DRHO_LEN_SPACE - 1);
            RHS_rho(sequence)
                << invertDiagonal(app.computeQ(i)) * app.GrhoL(sequence);
        }
        app.RHS = -(D1 * RHS_m + D2 * RHS_rho - app.GlambdaL);
        // At last, run the parallel-in-time TriMGRIT simulation
        core.Drive();

        // Collect resultant vectors for dlambda from the cores across which the
        // simulation was distributed, compute line search, and push results
        // back to parallel cores.
        if (rank == 0) {
            // Main processor
            // Receive from access the results of workers' computation
            std::vector<Vector> dlambda(DLAMBDA_LEN_TIME);

            int num_received = 0;
            int num_messages_expected = DLAMBDA_LEN_TIME;

            while (num_received < num_messages_expected) {
                int index;
                int message_size =
                    sizeof(int) + sizeof(double) * DLAMBDA_LEN_SPACE;
                char *buffer = (char *)malloc(message_size);
                MPI_Recv((void *)buffer, message_size, MPI_BYTE, MPI_ANY_SOURCE,
                         TAG_WORKER_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // printf("received: %d out of %d\n", num_received,
                // num_messages_expected);
                int *buffer_ = (int *)buffer;
                index = *(buffer_++);
                double *dbuffer = (double *)buffer_;
                Vector dlambda_val(DLAMBDA_LEN_SPACE);
                for (int i = 0; i < DLAMBDA_LEN_SPACE; i++) {
                    dlambda_val[i] = dbuffer[i];
                }
                dlambda[index] = dlambda_val;
                num_received++;
                free(buffer);
            }
            // Compute global new m, rho, and lambda using line search etc.
            Vector dlambda_long(dlambda.size() * dlambda[0].size());
            for (unsigned long i = 0; i < dlambda.size(); i++) {
                for (int j = 0; j < dlambda[0].size(); j++) {
                    dlambda_long[i * dlambda[0].size() + j] = dlambda[i][j];
                }
            }
            Vector dm_RHS = D1.transpose() * dlambda_long + app.GmL;
            Vector dm_long(DM_LEN_TIME * DM_LEN_SPACE);
            for (int i = 0; i < DM_LEN_TIME; i++) {
                dm_long(
                    Eigen::seq(i * DM_LEN_SPACE, (i + 1) * DM_LEN_SPACE - 1))
                    << -1.0 * invertDiagonal(app.computeP(i)) *
                           dm_RHS(Eigen::seq(i * DM_LEN_SPACE,
                                             (i + 1) * DM_LEN_SPACE - 1));
            }
            Vector drho_RHS = D2.transpose() * dlambda_long + app.GrhoL;
            Vector drho_long(DRHO_LEN_TIME * DRHO_LEN_SPACE);
            for (int i = 0; i < DRHO_LEN_TIME; i++) {
                drho_long(Eigen::seq(i * DRHO_LEN_SPACE,
                                     (i + 1) * DRHO_LEN_SPACE - 1))
                    << -1.0 * invertDiagonal(app.computeQ(i)) *
                           drho_RHS(Eigen::seq(i * DRHO_LEN_SPACE,
                                               (i + 1) * DRHO_LEN_SPACE - 1));
            }

            double alpha =
                line_search(dm_long, drho_long, dlambda_long, m_long, rho_long,
                            lambda_long, As, At, D1, D2, q_long, D);

            printf("Line search completed for iteration %d. alpha = %f\n", i_,
                   alpha);

            for (unsigned long i = 0; i < app.m.size(); i++) {
                app.m[i] +=
                    alpha * dm_long(Eigen::seq(i * DM_LEN_SPACE,
                                               (i + 1) * DM_LEN_SPACE - 1));
            }
            for (unsigned long i = 0; i < app.rho.size(); i++) {
                app.rho[i] +=
                    alpha * drho_long(Eigen::seq(i * DRHO_LEN_SPACE,
                                                 (i + 1) * DRHO_LEN_SPACE - 1));
            }
            for (unsigned long i = 0; i < app.lambda.size(); i++) {
                app.lambda[i] += alpha * dlambda[i];
            }
            // Send global picture of m, rho, and lambda
            int message_length =
                sizeof(double) *
                (DM_LEN_SPACE * DM_LEN_TIME + DRHO_LEN_SPACE * DRHO_LEN_TIME +
                 DLAMBDA_LEN_SPACE * DLAMBDA_LEN_TIME);
            double *buffer = (double *)malloc(message_length * sizeof(double));
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
            // Worker access should have already sent their results to the
            // main processor
            // Receive a global picture of m, rho, and lambda
            int message_length =
                sizeof(double) *
                (DM_LEN_SPACE * DM_LEN_TIME + DRHO_LEN_SPACE * DRHO_LEN_TIME +
                 DLAMBDA_LEN_SPACE * DLAMBDA_LEN_TIME);
            double *buffer = (double *)malloc(message_length * sizeof(double));
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

        // Decrease the normalization coefficient by a factor of 2.
        app.normcoeff /= 2;
    }

    for (unsigned long i = 0; i < app.rho.size(); i++) {
        std::cout << app.rho[i] << "\n" << std::endl;
    }

    // Print runtime to file (for runtime comparisons)
    // time = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Total Run Time: %f s \n", time);

    MPI_Finalize();

    return 0;
}
