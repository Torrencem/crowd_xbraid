/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory. Written by
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
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
 * Example:       ex-01.c
 *
 * Interface:     C
 *
 * Requires:      only C-language support
 *
 * Compile with:  make ex-01
 *
 * Help with:     this is the simplest example available, read the source
 *
 * Sample run:    mpirun -np 2 ex-01
 *
 * Description:   solve the scalar ODE
 *                   u' = lambda u,
 *                   with lambda=-1 and y(0) = 1
 *                in a very simplified XBraid setting.
 *
 *                When run with the default 10 time steps, the solution is:
 *                $ ./ex-01
 *                $ cat ex-01.out.00*
 *                  1.00000000000000e+00
 *                  6.66666666666667e-01
 *                  4.44444444444444e-01
 *                  2.96296296296296e-01
 *                  1.97530864197531e-01
 *                  1.31687242798354e-01
 *                  8.77914951989026e-02
 *                  5.85276634659351e-02
 *                  3.90184423106234e-02
 *                  2.60122948737489e-02
 *                  1.73415299158326e-02
 **/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "_braid.h"

#define TAG_DATA_TO_MASTER 2
#define TAG_LINE_SEARCH_RESULT 3
#define TAG_US_PREV 4

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct {
    int rank;
    int first_access;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct {
    double value;
} my_Vector;

int my_Step(braid_App app, braid_Vector ustop, braid_Vector fstop,
            braid_Vector u, braid_StepStatus status) {
    double tstart; /* current time */
    double tstop;  /* evolve to this time*/
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

    /* Use backward Euler to propagate solution */
    (u->value) = 1. / (1. + tstop - tstart) * (u->value);

    return 0;
}

int my_Init(braid_App app, double t, braid_Vector *u_ptr) {
    my_Vector *u;

    u = (my_Vector *)malloc(sizeof(my_Vector));
    if (t == 0.0) /* Initial condition */
    {
        (u->value) = 1.0;
    } else /* All other time points set to arbitrary value */
    {
        (u->value) = 0.456;
    }
    *u_ptr = u;

    return 0;
}

int my_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr) {
    my_Vector *v;

    v = (my_Vector *)malloc(sizeof(my_Vector));
    (v->value) = (u->value);
    *v_ptr = v;

    return 0;
}

int my_Free(braid_App app, braid_Vector u) {
    free(u);
    return 0;
}

int my_Sum(braid_App app, double alpha, braid_Vector x, double beta,
           braid_Vector y) {
    (y->value) = alpha * (x->value) + beta * (y->value);
    return 0;
}

int my_SpatialNorm(braid_App app, braid_Vector u, double *norm_ptr) {
    double dot;

    dot = (u->value) * (u->value);
    *norm_ptr = sqrt(dot);
    return 0;
}

// --- LINE SEARCH STUFF ---

double objective(double alpha, braid_BaseVector *us_prev, braid_BaseVector *us, int len, int level,
                 braid_App app, braid_Core core) {
    int cfactor;
    _braid_GetCFactor(core, level, &cfactor);
    double result = 0.0;
    // r is reused on each iteration. It stores the residual calculation
    struct _braid_BaseVector_struct *r;
    _braid_BaseInit(core, app, 0.0, &r);
    for (int i = 0; i < len; i++) {
        struct _braid_BaseVector_struct *us_updated;
        _braid_BaseClone(core, app, us[i], &us_updated);
        _braid_BaseSum(core, app, 1 - alpha, us_prev[i], alpha, us_updated);
        // Level + 1 since we're doing the line search on the coarse grid
        // Compute residual manually as ustop - Phi(ustart)
        /* _braid_Residual(core, level, i, ustop, r); */
        {
           braid_Real       tol      = _braid_CoreElt(core, tol);
           braid_Int        iter     = _braid_CoreElt(core, niter);
           _braid_Grid    **grids    = _braid_CoreElt(core, grids);
           braid_StepStatus status   = (braid_StepStatus)core;
           braid_Int        nrefine  = _braid_CoreElt(core, nrefine);
           braid_Int        gupper   = _braid_CoreElt(core, gupper);
           braid_Int        ilower   = _braid_GridElt(grids[level], ilower);
           braid_Real      *ta       = _braid_GridElt(grids[level], ta);

           braid_Int        ii;

           ii = i-ilower;
           _braid_StepStatusInit(ta[ii-1], ta[ii], i-1, tol, iter, level, nrefine, gupper, status);

           // Now compute Phi(us_prev[i])
        
           braid_BaseVector phi_u_prev;
           _braid_BaseClone(core, app, us_prev[i], &phi_u_prev);
            
           // Workaround for a bug in xbraid
            struct _braid_BaseVector_struct ustop;
            ustop.userVector = NULL;
            ustop.bar = NULL;

            _braid_BaseStep(core, app, &ustop, NULL, phi_u_prev, level + 1, status);
            // r = ustop - \Phi(ustart)
            _braid_BaseSum(core, app, 1.0, us[i], -1.0, phi_u_prev);
            r = phi_u_prev;
        }
        double norm;
        _braid_BaseSpatialNorm(core, app, r, &norm);
        result += norm;
        _braid_BaseFree(core, app, us_updated);
    }
    _braid_BaseFree(core, app, r);
    return result;
}

typedef double (*Objective)(double, braid_BaseVector *, braid_BaseVector *, int, int, braid_App, braid_Core);

double line_search(Objective f, braid_BaseVector *us_prev, braid_BaseVector *us, int len, int level,
                   braid_App app, braid_Core status) {
    /* return 0.0; */
    double a = 0.0;
    double b = 1.0;
    double gr = (1.0 + sqrt(5.0)) / 2.0;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;
    while (fabs(c - d) > 0.00001) {
        if (f(c, us_prev, us, len, level, app, status) < f(d, us_prev, us, len, level, app, status)) {
            b = d;
        } else {
            a = c;
        }
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    }
    return (a + b) / 2;
}

// --- END LINE SEARCH STUFF ---

int my_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus) {
    return 0;
}

int my_Sync(braid_App app, braid_SyncStatus status) {
    MPI_Status mpi_status;
    MPI_Request request;
    int num_vectors;
    int level;
    // The lower and upper indices in time that our proccess has computed
    int lower_t, upper_t;
    int num_procesors;
    int num_total_points, c_factor;
    braid_Core core = (braid_Core) status;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procesors);

    braid_SyncStatusGetLevel(status, &level);

    braid_SyncStatusGetNTPoints(status, &num_total_points);
    _braid_GetCFactor(core, level, &c_factor);
    
    if (level != 0) {
        return 0;
    }

    num_vectors = num_total_points / c_factor;
    
    // Collect us - the vectors for our process
    _braid_Grid **grids = _braid_StatusElt(status, grids);
    _braid_Grid *grid = grids[level];
    upper_t = grid->cupper / c_factor;
    lower_t = grid->clower / c_factor;

    int bvector_size; 
    _braid_BaseBufSize(core, app, &bvector_size, (braid_BufferStatus) core);
    int safe_message_size =
        (num_vectors / num_procesors + 2) * bvector_size + 2 * sizeof(int);

    int my_num_values = upper_t - lower_t; // Assuming upper_t is not inclusive

    assert(2 * sizeof(int) + my_num_values * bvector_size <= (unsigned long) safe_message_size);
    
    braid_BaseVector *us = malloc(my_num_values * my_num_values);

    for (int i = 0; i < my_num_values; i++) {
        int fine_index;
        _braid_MapCoarseToFine(lower_t + i, c_factor, fine_index);
        braid_BaseVector base_v;

        _braid_UGetVectorRef((braid_Core) status, level, fine_index, &base_v);

        assert(base_v != NULL);

        us[i] = base_v;
    }

    // Now that we have us...

    if (app->rank == 0) {
        // Master process
        // Receive all other vectors
        braid_BaseVector *us_combined = malloc(sizeof(braid_BaseVector) * num_vectors);
        // handle our own us
        for (int i = lower_t; i < upper_t; i++) {
            us_combined[lower_t + i] = us[i];
        }
        for (int i = 0; i < num_procesors - 1; i++) {
            // Received message:
            // (int: the start index (lower_t)) : (int: the number of us (upper_t - lower_t)) : values (us)
            char *message = malloc(safe_message_size);
            MPI_Recv(message, safe_message_size, MPI_BYTE, MPI_ANY_SOURCE, TAG_DATA_TO_MASTER, MPI_COMM_WORLD, &mpi_status);
            // Finagle the message
            int start_index = ((int *) message)[0];
            int us_len = ((int *) message)[1];
            char *dmessage = (char *) (message + 2 * sizeof(int));
            for (int i = 0; i < us_len; i++) {
                /* us_combined[i + start_index] = dmessage[i]; */
                // use BufUnpack
                int index = i + start_index;
                _braid_BaseBufUnpack(core, app, dmessage + i * bvector_size, us_combined + index, (braid_BufferStatus) core);
            }
            free(dmessage);
        }
        // Perform the line search
        if (app->first_access) {
            app->first_access = 0;
        } else {
            char *us_prev_message = malloc(bvector_size * num_vectors);
            // Receive us_prev from previous iteration of access
            MPI_Recv(us_prev_message, num_vectors * bvector_size, MPI_BYTE, 0, TAG_US_PREV,
                     MPI_COMM_WORLD, &mpi_status);
            // Unpack from message into us_prev
            braid_BaseVector *us_prev = malloc(sizeof(braid_BaseVector) * num_vectors);
            for (int i = 0; i < num_vectors; i++) {
                _braid_BaseBufUnpack(core, app, us_prev_message + i * bvector_size, us_prev + i, (braid_BufferStatus) core);
            }

            free(us_prev_message);

            double alpha =
                line_search(objective, us_prev, us_combined, num_vectors, level, app, (braid_Core) status);
            // ... put result in us_combined
            for (int i = 0; i < num_vectors; i++) {
                _braid_BaseSum(core, app, 1 - alpha, us_prev[i], alpha, us_combined[i]);
            }

            // As the main process, we need to update our us[i]'s
            for (int i = lower_t; i < upper_t; i++) {
                _braid_BaseClone(core, app, us_combined[i], us + i - lower_t);
            }
        }
        // Send us_prev for next iteration of loop
        // Pack us_combined to send
        char *us_combined_packed = malloc(num_vectors * bvector_size);
        for (int i = 0; i < num_vectors; i++) {
            _braid_BaseBufPack(core, app, us_combined[i], us_combined_packed + i * bvector_size, (braid_BufferStatus) core);
        }
        MPI_Isend(us_combined_packed, num_vectors * bvector_size, MPI_BYTE, 0, TAG_US_PREV, MPI_COMM_WORLD, &request);
        // Send updated us_combined to workers
        MPI_Bcast(us_combined_packed, num_vectors * bvector_size, MPI_BYTE, 0, MPI_COMM_WORLD);
        free(us_combined);
        free(us_combined_packed);
    } else {
        // Worker process
        // Send u's data to master process
        char *message = malloc(safe_message_size);
        ((int *) message)[0] = lower_t;
        ((int *) message)[1] = my_num_values;
        char *dmessage = (char *) (message + 2 * sizeof(int));
        for (int i = 0; i < my_num_values; i++) {
            /* dmessage[i] = us[i]->value; */
            // Use BufPack
            _braid_BaseBufPack(core, app, us[i], dmessage + i * bvector_size, (braid_BufferStatus) core);
        }
        MPI_Send(message, safe_message_size, MPI_BYTE, 0, TAG_DATA_TO_MASTER, MPI_COMM_WORLD);
        free(message);
        // Receive updated us values
        char *updated_us = malloc(num_vectors * bvector_size);
        MPI_Bcast(updated_us, num_vectors, MPI_BYTE, 0, MPI_COMM_WORLD);
        // Put them back into us
        for (int i = 0; i < my_num_values; i++) {
            /* us[i]->value = updated_us[i + lower_t]; */
            // Use BufUnpack
            _braid_BaseBufUnpack(core, app, updated_us + i * bvector_size, us + i, (braid_BufferStatus) core);
        }
        free(updated_us);
    }

    free(us);

    return 0;
}

int my_BufSize(braid_App app, int *size_ptr, braid_BufferStatus bstatus) {
    *size_ptr = sizeof(double);
    return 0;
}

int my_BufPack(braid_App app, braid_Vector u, void *buffer,
               braid_BufferStatus bstatus) {
    double *dbuffer = buffer;

    dbuffer[0] = (u->value);
    braid_BufferStatusSetSize(bstatus, sizeof(double));

    return 0;
}

int my_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr,
                 braid_BufferStatus bstatus) {
    double *dbuffer = buffer;
    my_Vector *u;

    u = (my_Vector *)malloc(sizeof(my_Vector));
    (u->value) = dbuffer[0];
    *u_ptr = u;

    return 0;
}

double my_SpaceTimeNorm(braid_App app, int time_steps, braid_Vector *u_ptr) {
    double accumulator = 0.0;
    double norm;
    for (int i = 0; i < time_steps; i++) {
        my_SpatialNorm(app, u_ptr[i], &norm);
        accumulator += (norm * norm);
    }
    return sqrt(norm);
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    braid_Core core;
    my_App *app;
    double tstart, tstop;
    int ntime, rank;

    /* Define time domain: ntime intervals */
    ntime = 10;
    tstart = 0.0;
    tstop = tstart + ntime / 2.;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* set up app structure */
    app = (my_App *)malloc(sizeof(my_App));
    (app->rank) = rank;
    app->first_access = 1;

    /* initialize XBraid and set options */
    braid_Init(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
               my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm,
               my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

    /* Set some typical Braid parameters */
    braid_SetPrintLevel(core, 2);
    braid_SetMaxLevels(core, 4);
    braid_SetAbsTol(core, 1.0e-06);
    braid_SetCFactor(core, -1, 2);

    /* For our custom line searching */
    braid_SetSync(core, my_Sync);

    /* Run simulation, and then clean up */
    braid_Drive(core);

    braid_Destroy(core);
    free(app);
    MPI_Finalize();

    return (0);
}
