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

// Crowd model problem
//
// This program simulates an optimization problem of modelling the movement of a
// crowd from one configuration to another. The math behind the code used here
// can be found in Haber and Horesh's 2015 paper "A Multilevel Method for the
// Solution of Time Dependent Optimal Transport".  We use the same algorithm
// as in crowd_horesh.cpp, but include TriMGRIT as our matrix equation solver.

#define CROWD_HORESH_LIBRARY
#include "crowd_horesh.cpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "braid_test.h"
#include "tribraid.hpp"
#define PI 3.14159265
#define TAG_WORKER_RESULT 42
#define TYPE_RHO 2
#define TYPE_LAMBDA 3
#define TYPE_M 4

#define DM_LEN_SPACE (mspace + 1)
#define DRHO_LEN_SPACE (mspace)
#define DLAMBDA_LEN_SPACE (mspace)
#define Q_LEN_SPACE (mspace)

#define DM_LEN_TIME (ntime)
#define DRHO_LEN_TIME (ntime + 1)
#define DLAMBDA_LEN_TIME (ntime + 2)
#define Q_LEN_TIME (ntime + 2)

Sparse invertDiagonal(const Sparse &a) {
    // assert(a.rows() == a.cols());
    Sparse m2(a.rows(), a.cols());

    for (int i = 0; i < a.rows(); i++) {
        m2.insert(i, i) = 1.0 / a.coeff(i, i);
    }

    return m2;
}

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

class BraidVector {
  public:
    int index = -1;
    Vector dm, drho, dlambda;

    BraidVector(Vector dm_, Vector drho_, Vector dlambda_)
        : dm(dm_), drho(drho_), dlambda(dlambda_) {}

    BraidVector(int mspace) {
        this->dm = Vector(DM_LEN_SPACE);
        this->dm.setZero();

        this->drho = Vector(DRHO_LEN_SPACE);
        this->drho.setZero();

        this->dlambda = Vector(DLAMBDA_LEN_SPACE);
        this->dlambda.setZero();
    }

    virtual ~BraidVector(){};
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

    // Not needed
    virtual int Residual(braid_Vector u_, braid_Vector r_,
                         BraidStepStatus &pstatus) override {
        return 0;
    }

    const Vector apply_Phi(Vector &u, Vector &v, const double t,
                           const double dt);

    // Compute \nabla_m_i
    const Vector compute_GwMi(int index);
    // Compute \nabla_\lambda_i
    const Vector compute_GwLi(int index);
    // Compute \nabla_\rho_{i - 1/2}
    const Vector compute_GwRhoi(int index);
};

const Vector MyBraidApp::compute_GwMi(int index) {
    Eigen::Map<Vector> res(this->GmL.data() + (index * DM_LEN_SPACE),
                           DM_LEN_SPACE);

    Vector res_(res);
    return res_;
}

const Vector MyBraidApp::compute_GwLi(int index) {
    Eigen::Map<Vector> res(this->GlambdaL.data() + (index * DLAMBDA_LEN_SPACE),
                           DLAMBDA_LEN_SPACE);

    Vector res_(res);
    return res_;
}

const Vector MyBraidApp::compute_GwRhoi(int index) {
    Eigen::Map<Vector> res(this->GrhoL.data() + (index * DRHO_LEN_SPACE),
                           DRHO_LEN_SPACE);

    Vector res_(res);
    return res_;
}

MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_,
                       double tstop_, int ntime_)
    : TriBraidApp(comm_t_, tstart_, tstop_, ntime_) {}

/*------------------------------------*/

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute r = A(r) - f */

int MyBraidApp::TriResidual(braid_Vector uleft_, braid_Vector uright_,
                            braid_Vector f_, braid_Vector r_,
                            BraidTriStatus &status) {
    BraidVector *uleft = (BraidVector *)uleft_;
    BraidVector *uright = (BraidVector *)uright_;
    BraidVector *f = (BraidVector *)f_;
    BraidVector *r = (BraidVector *)r_;

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

    Sparse X((mspace + 1), mspace);

    for (int i = 0; i < mspace; i++) {
        X.insert(i, i) = 0.25;
        X.insert(i + 1, i) = 0.25;
    }

    Vector mi;

    if (index != 0 && index != DM_LEN_TIME + 1) {
        mi = m[index];
    } else {
        mi = Vector(DM_LEN_SPACE);
        mi.setZero();
    }

    Vector mim1(DM_LEN_SPACE);

    if (index == 0 || index == 1) {
        mim1.setZero();
    } else {
        mim1 = m[index - 1];
    }

    // Compute Qi and Pi
    // Vector tmp1 = rhoi.cwiseInverse() + rhoip1.cwiseInverse();
    Vector tmp1(DRHO_LEN_SPACE);
    tmp1.setZero();
    if (index != 0) {
        // Able to add rhoi
        tmp1 += rho[index].cwiseInverse();
    }
    if (index != DRHO_LEN_TIME) { // TODO: Is this condition correct?
        // Able to add rhoip1
        tmp1 += rho[index + 1].cwiseInverse();
    }
    Vector tmp1_ = X * tmp1;
    auto Pi_ = 2.0 * (tmp1_).asDiagonal();
    assert(mi.cwiseProduct(mi).size() == mim1.cwiseProduct(mim1).size());
    Vector tmp2 = mi.cwiseProduct(mi) + mim1.cwiseProduct(mim1);
    Vector tmp3;
    if (index != 0 && index < DRHO_LEN_TIME) {
        Vector rhoi = rho[index];
        tmp3 = rhoi.cwiseProduct(rhoi.cwiseProduct(rhoi)).cwiseInverse();
    } else {
        tmp3 = Vector(DRHO_LEN_SPACE);
        tmp3.setZero();
    }
    Vector tmp4 = X.transpose() * tmp2;
    Vector Qi_(tmp3.size());
    Qi_.setConstant(2.0);
    Qi_ = Qi_.cwiseProduct(tmp3).cwiseProduct(tmp4);
    auto Qi__ = Qi_.asDiagonal();
    Sparse Qi(Qi__);
    Sparse Pi(Pi_);

    Vector nabla_m_i = this->compute_GwMi(index);
    Vector nabla_rho_i = this->compute_GwRhoi(index);
    Vector nabla_lambda_i = this->compute_GwLi(index);

    Vector dlambda_left(DLAMBDA_LEN_SPACE);
    if (uleft != nullptr) {
        dlambda_left = uleft->dlambda;
    } else {
        dlambda_left.setZero();
    }

    Vector drho_right(DRHO_LEN_SPACE);
    if (uright != nullptr) {
        drho_right = uright->drho;
    } else {
        drho_right.setZero();
    }

    /*if (f == nullptr) {
        BraidVector f(mspace);
        r->dm = Pi * f.dm + K.transpose() * f.dlambda + nabla_m_i;
        Vector tmp1 = Qi * f.drho + dlambda_left / dt - f.dlambda / dt;
        r->drho = tmp1 + nabla_rho_i;
        r->dlambda = K * f.dm + drho_right / dt - f.drho / dt + nabla_lambda_i;
    } else {
        r->dm = Pi * f->dm + K.transpose() * f->dlambda + nabla_m_i;
        r->drho =
            Qi * f->drho + dlambda_left / dt - f->dlambda / dt + nabla_rho_i;
        r->dlambda =
            K * f->dm + drho_right / dt - f->drho / dt + nabla_lambda_i;
    }*/
    r->dm = Pi * r->dm + K.transpose() * r->dlambda;
    r->drho = Qi * r->drho + dlambda_left / dt - r->dlambda / dt;
    r->dlambda = K * r->dm + drho_right / dt - r->drho / dt;
    if (f != nullptr) {
        r->dm = r->dm - f->dm;
        r->drho = r->drho - f->drho;
        r->dlambda = r->dlambda - f->dlambda;
    }
    return 0;
}

/*------------------------------------*/

/* Solve A(u) = f */

int MyBraidApp::TriSolve(braid_Vector uleft_, braid_Vector uright_,
                         braid_Vector fleft_, braid_Vector fright_,
                         braid_Vector f_, braid_Vector u_,
                         BraidTriStatus &status) {
    BraidVector *uleft = (BraidVector *)uleft_;
    BraidVector *uright = (BraidVector *)uright_;
    // BraidVector *f = (BraidVector *) f_;
    BraidVector *u = (BraidVector *)u_;
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
    //    final_index -= 1;

    //    std::cout << "Index " << index << " and final index " << final_index
    //    << std::endl;
    Sparse X((mspace + 1), mspace);

    for (int i = 0; i < mspace; i++) {
        X.insert(i, i) = 0.25;
        X.insert(i + 1, i) = 0.25;
    }

    Vector mi;

    if (index != 0 && index != DM_LEN_TIME + 1) {
        mi = m[index];
    } else {
        mi = Vector(DM_LEN_SPACE);
        mi.setZero();
    }

    Vector mim1(DM_LEN_SPACE);

    if (index == 0 || index == 1) {
        mim1.setZero();
    } else {
        mim1 = m[index - 1];
    }

    // Compute Qi and Pi
    // Vector tmp1 = rhoi.cwiseInverse() + rhoip1.cwiseInverse();
    Vector tmp1(DRHO_LEN_SPACE);
    tmp1.setZero();
    if (index != 0) {
        // Able to add rhoi
        tmp1 += rho[index].cwiseInverse();
    }
    if (index != DRHO_LEN_TIME) {
        // Able to add rhoip1
        tmp1 += rho[index + 1].cwiseInverse();
    }
    Vector tmp1_ = X * tmp1;
    auto Pi_ = 2.0 * (tmp1_).asDiagonal();
    assert(mi.cwiseProduct(mi).size() == mim1.cwiseProduct(mim1).size());
    Vector tmp2 = mi.cwiseProduct(mi) + mim1.cwiseProduct(mim1);
    Vector tmp3;
    if (index != 0 && index < DRHO_LEN_TIME) {
        Vector rhoi = rho[index];
        tmp3 = rhoi.cwiseProduct(rhoi.cwiseProduct(rhoi)).cwiseInverse();
    } else {
        tmp3 = Vector(DRHO_LEN_SPACE);
        tmp3.setZero();
    }
    Vector tmp4 = X.transpose() * tmp2;
    Vector Qi_(tmp3.size());
    Qi_.setConstant(2.0);
    Qi_ = Qi_.cwiseProduct(tmp3).cwiseProduct(tmp4);
    auto Qi__ = Qi_.asDiagonal();
    Sparse Qi(Qi__);
    Sparse Pi(Pi_);

    if (index == 1) {
        // Compute delta_rho_1
        // Equation 12
        Vector nabla_l_0 = this->compute_GwLi(0);
        u->drho = dt * -nabla_l_0;

        // Compute delta_lambda_1
        // Equation 13
        Vector nabla_rho_1 = this->compute_GwRhoi(0);
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
        if (index != 0) {
            Vector nabla_m_i = this->compute_GwMi(index);
            Vector nabla_lambda_i = this->compute_GwLi(index);
            Vector nabla_rho_i = this->compute_GwRhoi(index);
            Vector delta_lambda_im1;
            delta_lambda_im1 = uleft->dlambda;
            Vector delta_rho_ip1;
            if (uright == nullptr) {
                delta_rho_ip1 = Vector(u->drho.size());
                delta_rho_ip1.setZero();
            } else {
                delta_rho_ip1 = uright->drho;
            }
            Sparse I(DLAMBDA_LEN_SPACE, DLAMBDA_LEN_SPACE);
            I.setIdentity();
            Sparse A = -dt * Qi * K * invertDiagonal(Pi) * K.transpose();
            A -= I / dt;
            Vector b = dt * Qi * K * invertDiagonal(Pi) * nabla_m_i -
                       Qi * delta_rho_ip1 - dt * Qi * nabla_lambda_i -
                       delta_lambda_im1 / dt - nabla_rho_i;
            Eigen::BiCGSTAB<Sparse, Eigen::IncompleteLUT<double>> solver;
            solver.compute(A);
            u->dlambda = solver.solve(b);
            if (index != final_index) { // Otherwise delta_m doesn't matter
                Vector nabla_m_i = this->compute_GwMi(index);
                // Setup Ax = b system
                Vector b = -nabla_m_i - K.transpose() * u->dlambda;

                // Solve Pi x = b
                u->dm = invertDiagonal(Pi) * b;

                u->drho = dt * (nabla_lambda_i + K * u->dm +
                                (1 / dt) * delta_rho_ip1);
            } else {
                u->drho = dt * nabla_lambda_i;
            }

        } else {
            Vector nabla_lambda_0 = this->compute_GwLi(0);
            Vector dlambda_1 = uright->dlambda;
            Vector nabla_rho_1 = this->compute_GwRhoi(0);
            u->dlambda =
                dt * (dt * Qi * nabla_lambda_0 + dlambda_1 / dt - nabla_rho_1);
        }
        // Compute delta_m
        // Equation 9
    }
    /* no refinement */
    status.SetRFactor(1);

    return 0;
}

/*------------------------------------*/

/* This is only called from level 0 */

int MyBraidApp::Init(double t, braid_Vector *u_ptr_) {
    BraidVector **u_ptr = (BraidVector **)u_ptr_;

    Vector dm(DM_LEN_SPACE);
    Vector drho(DRHO_LEN_SPACE);
    Vector dlambda(DLAMBDA_LEN_SPACE);

    dm.setConstant(0.0);
    drho.setConstant(0.0);
    dlambda.setConstant(0.0);

    dlambda.setRandom();

    *u_ptr = new BraidVector(dm, drho, dlambda);

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Clone(braid_Vector u_, braid_Vector *v_ptr_) {
    BraidVector *u = (BraidVector *)u_;
    BraidVector **v_ptr = (BraidVector **)v_ptr_;

    // TODO: maybe?:
    // delete *v_ptr;

    *v_ptr = new BraidVector(u->dm, u->drho, u->dlambda);
    (*v_ptr)->index = u->index;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Free(braid_Vector u_) {
    BraidVector *u = (BraidVector *)u_;
    delete u;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Sum(double alpha, braid_Vector x_, double beta,
                    braid_Vector y_) {
    BraidVector *y = (BraidVector *)y_;
    BraidVector *x = (BraidVector *)x_;

    y->dm = alpha * x->dm + beta * y->dm;
    y->drho = alpha * x->drho + beta * y->drho;
    y->dlambda = alpha * x->dlambda + beta * y->dlambda;

    return 0;
}

/*------------------------------------*/

int MyBraidApp::SpatialNorm(braid_Vector u_, double *norm_ptr) {
    BraidVector *u = (BraidVector *)u_;

    *norm_ptr = sqrt(u->dlambda.squaredNorm() + u->dm.squaredNorm() +
                     u->drho.squaredNorm());

    return 0;
}

/*------------------------------------*/

int MyBraidApp::Access(braid_Vector u_, BraidAccessStatus &astatus) {
    BraidVector *u = (BraidVector *)u_;

    int done, index;

    astatus.GetDone(&done);
    if (done) {
        astatus.GetTIndex(&index);
        std::cout << "At index " << index << std::endl;
        std::cout << "drho is " << u->drho << std::endl;
        std::cout << "dm is " << u->dm << std::endl;
        std::cout << "dlambda is " << u->dlambda << std::endl;

        MPI_Request send_request;
        int message_size = sizeof(int) * 2 + sizeof(double) * (mspace + 2);

        // type is TYPE_RHO
        if (index != 0) {
            char *buffer = (char *)calloc(message_size, 1);
            int *ibuffer = (int *)buffer;
            *(ibuffer++) = TYPE_RHO;
            *(ibuffer++) = index;
            double *dbuffer = (double *)ibuffer;
            for (int i = 0; i < DRHO_LEN_SPACE; i++) {
                dbuffer[i] = u->drho[i];
            }
            MPI_Isend((void *)buffer, message_size, MPI_BYTE, 0,
                      TAG_WORKER_RESULT, MPI_COMM_WORLD, &send_request);
            free(buffer);
        }
        // type is TYPE_LAMBDA
        {
            char *buffer = (char *)calloc(message_size, 1);
            int *ibuffer = (int *)buffer;
            *(ibuffer++) = TYPE_LAMBDA;
            *(ibuffer++) = index;
            double *dbuffer = (double *)ibuffer;
            for (int i = 0; i < DLAMBDA_LEN_SPACE; i++) {
                dbuffer[i] = u->dlambda[i];
            }
            MPI_Isend((void *)buffer, message_size, MPI_BYTE, 0,
                      TAG_WORKER_RESULT, MPI_COMM_WORLD, &send_request);
            free(buffer);
        }
        // type is TYPE_M
        if (index != 0 && index != DM_LEN_TIME + 1) {
            char *buffer = (char *)calloc(message_size, 1);
            int *ibuffer = (int *)buffer;
            *(ibuffer++) = TYPE_M;
            *(ibuffer++) = index;
            double *dbuffer = (double *)ibuffer;
            for (int i = 0; i < DM_LEN_SPACE; i++) {
                dbuffer[i] = u->dm[i];
            }
            MPI_Isend((void *)buffer, message_size, MPI_BYTE, 0,
                      TAG_WORKER_RESULT, MPI_COMM_WORLD, &send_request);
            free(buffer);
        }
    }

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufSize(int *size_ptr, BraidBufferStatus &bstatus) {
    // sizeof(index) + sizeof(drho) + sizeof(dlambda) + sizeof(dm)
    *size_ptr = sizeof(int) + DRHO_LEN_SPACE * sizeof(double) +
                DLAMBDA_LEN_SPACE * sizeof(double) +
                DM_LEN_SPACE * sizeof(double);
    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufPack(braid_Vector u_, void *buffer_,
                        BraidBufferStatus &bstatus) {
    BraidVector *u = (BraidVector *)u_;
    int *buffer = (int *)buffer_;

    *(buffer++) = u->index;

    double *dbuffer = (double *)buffer;

    for (int i = 0; i < DRHO_LEN_SPACE; i++, dbuffer++) {
        *dbuffer = u->drho[i];
    }

    for (int i = 0; i < DLAMBDA_LEN_SPACE; i++, dbuffer++) {
        *dbuffer = u->dlambda[i];
    }

    for (int i = 0; i < DM_LEN_SPACE; i++, dbuffer++) {
        *dbuffer = u->dm[i];
    }

    bstatus.SetSize(sizeof(int) + DRHO_LEN_SPACE * sizeof(double) +
                    DLAMBDA_LEN_SPACE * sizeof(double) +
                    DM_LEN_SPACE * sizeof(double));

    return 0;
}

/*------------------------------------*/

int MyBraidApp::BufUnpack(void *buffer_, braid_Vector *u_ptr_,
                          BraidBufferStatus &bstatus) {
    BraidVector **u_ptr = (BraidVector **)u_ptr_;

    int *buffer = (int *)buffer_;

    (*u_ptr)->index = *(buffer++);

    double *dbuffer = (double *)buffer;

    for (int i = 0; i < DRHO_LEN_SPACE; i++, dbuffer++) {
        (*u_ptr)->drho[i] = *dbuffer;
    }

    for (int i = 0; i < DLAMBDA_LEN_SPACE; i++, dbuffer++) {
        (*u_ptr)->dlambda[i] = *dbuffer;
    }

    for (int i = 0; i < DM_LEN_SPACE; i++, dbuffer++) {
        (*u_ptr)->dm[i] = *dbuffer;
    }

    bstatus.SetSize(sizeof(int) + DRHO_LEN_SPACE * sizeof(double) +
                    DLAMBDA_LEN_SPACE * sizeof(double) +
                    DM_LEN_SPACE * sizeof(double));

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
    mspace = 4;
    ntime = 4;

    /* Define some Braid parameters */
    max_levels = 1;
    min_coarse = 1;
    nrelax = 10;
    nrelaxc = 10;
    maxiter = 5;
    cfactor = 2;
    tol = 1.0e-6;
    access_level = 1;
    print_level = 2;

    /* Define the space step */
    dx = (double)1 / (mspace + 1);

    /* Define time domain and step */
    tstart = 0.0; /* Beginning of time domain */
    tstop = 1.0;  /* End of time domain*/

    /* Set up the app structure */
    // ntime + 1 for the 2 extra time points caused by the staggered grid
    auto app = MyBraidApp(MPI_COMM_WORLD, rank, tstart, tstop, ntime + 1);
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
        core.SetNRelax(max_levels - 1, nrelaxc); /* nrelax on coarsest level */
    }
    core.SetCFactor(-1, cfactor);
    core.SetAccessLevel(access_level);
    core.SetPrintLevel(print_level);
    core.SetMaxIter(maxiter);
    core.SetAbsTol(tol);

    /* Parallel-in-time TriMGRIT simulation */
    time = 1.0;
    int iters = 1;

    double d_time = time / ntime;

    app.m = std::vector<Vector>();
    // M doesn't exist at the first time point
    app.m.push_back(Vector(DM_LEN_SPACE));
    for (int i = 0; i < DM_LEN_TIME; i++) {
        Vector m_val(DM_LEN_SPACE);
        m_val.setConstant(0.0);
        app.m.push_back(m_val);
    }

    app.rho = std::vector<Vector>();
    // Rho doesn't exist at the first time point
    app.rho.push_back(Vector(DRHO_LEN_SPACE));
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

    app.q = std::vector<Vector>(Q_LEN_TIME);

    double accumulator = 0.0;
    Vector q_val(Q_LEN_SPACE);
    for (int i = 0; i < mspace; i++) {
        q_val[i] = 0.5; // initial_condition(accumulator);
        accumulator += 1.0 / ((double)mspace - 1.0);
    }
    app.q[0] = q_val / d_time;

    for (int i = 1; i < ntime + 1; i++) {
        app.q[i] = q_val * 0.0;
    }

    accumulator = 0.0;
    q_val = Vector(mspace);
    for (int i = 0; i < mspace; i++) {
        q_val[i] = -0.5; // * final_condition(accumulator);
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

    Vector q_long(app.q.size() * app.q[0].size());
    for (unsigned long i = 0; i < app.q.size(); i++) {
        for (int j = 0; j < app.q[0].size(); j++) {
            double val = app.q[i][j];
            q_long[i * app.q[0].size() + j] = val;
        }
    }

    for (int i_ = 0; i_ < iters; i_++) {
        Vector rho_long((app.rho.size() - 1) * app.rho[1].size());
        for (unsigned long i = 1; i < app.rho.size(); i++) {
            for (int j = 0; j < app.rho[1].size(); j++) {
                rho_long[(i - 1) * app.rho[1].size() + j] = app.rho[i][j];
            }
        }
        Vector m_long((app.m.size() - 1) * app.m[1].size());
        for (unsigned long i = 1; i < app.m.size(); i++) {
            for (int j = 0; j < app.m[1].size(); j++) {
                m_long[(i - 1) * app.m[1].size() + j] = app.m[i][j];
            }
        }
        Vector lambda_long(app.lambda.size() * app.lambda[0].size());
        for (unsigned long i = 0; i < app.lambda.size(); i++) {
            for (int j = 0; j < app.lambda[0].size(); j++) {
                lambda_long[i * app.lambda[0].size() + j] = app.lambda[i][j];
            }
        }

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
                    D2.transpose() * lambda_long;

        core.Drive();

        if (rank == 0) {
            // Main processor
            // Receive from access the results of workers' computation
            std::vector<Vector> drho(DRHO_LEN_TIME);
            std::vector<Vector> dlambda(DLAMBDA_LEN_TIME);
            std::vector<Vector> dm(DM_LEN_TIME);

            int num_received = 0;
            int num_messages_expected =
                DRHO_LEN_TIME + DLAMBDA_LEN_TIME + DM_LEN_TIME;

            while (num_received < num_messages_expected) {
                int type, index;
                int message_size =
                    sizeof(int) * 2 + sizeof(double) * (mspace + 2);
                char *buffer = (char *)malloc(message_size);
                MPI_Recv((void *)buffer, message_size, MPI_BYTE, MPI_ANY_SOURCE,
                         TAG_WORKER_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // printf("received: %d out of %d\n", num_received,
                // num_messages_expected);
                int *buffer_ = (int *)buffer;
                type = *(buffer_++);
                index = *(buffer_++);
                double *dbuffer = (double *)buffer_;
                if (type == TYPE_RHO) {
                    Vector drho_val(DRHO_LEN_SPACE);
                    for (int i = 0; i < DRHO_LEN_SPACE; i++) {
                        drho_val[i] = dbuffer[i];
                    }
                    drho[index - 1] = drho_val;
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
                    dm[index - 1] = dm_val;
                }
                num_received++;
                free(buffer);
            }
            // Compute global new m, rho, and lambda using line search etc.
            Vector drho_long(drho.size() * drho[0].size());
            for (unsigned long i = 0; i < drho.size(); i++) {
                for (int j = 0; j < drho[0].size(); j++) {
                    drho_long[i * drho[0].size() + j] = drho[i][j];
                }
            }
            Vector dm_long(dm.size() * dm[0].size());
            for (unsigned long i = 0; i < dm.size(); i++) {
                for (int j = 0; j < dm[0].size(); j++) {
                    dm_long[i * dm[0].size() + j] = dm[i][j];
                }
            }
            Vector dlambda_long(dlambda.size() * dlambda[0].size());
            for (unsigned long i = 0; i < dlambda.size(); i++) {
                for (int j = 0; j < dlambda[0].size(); j++) {
                    dlambda_long[i * dlambda[0].size() + j] = dlambda[i][j];
                }
            }
            double alpha =
                line_search(dm_long, drho_long, dlambda_long, m_long, rho_long,
                            lambda_long, As, At, D1, D2, q_long, D);

            printf("Line search completed for iteration %d. alpha = %f\n", i_,
                   alpha);

            for (unsigned long i = 1; i < app.m.size(); i++) {
                app.m[i] += alpha * dm[i - 1];
            }
            for (unsigned long i = 1; i < app.rho.size(); i++) {
                app.rho[i] += alpha * drho[i - 1];
            }
            for (unsigned long i = 0; i < app.lambda.size(); i++) {
                app.lambda[i] += alpha * dlambda[i];
            }
            // Send global picture of m, rho, and lambda
            int message_length =
                sizeof(double) * (app.m.size() * app.m[0].size() +
                                  app.rho.size() * app.rho[0].size() +
                                  app.lambda.size() * app.lambda[0].size());
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
            // Worker
            // access should have already sent their results to the main
            // processor Receive a global picture of m, rho, and lambda
            int message_length =
                sizeof(double) * (app.m.size() * app.m[0].size() +
                                  app.rho.size() * app.rho[0].size() +
                                  app.lambda.size() * app.lambda[0].size());
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
