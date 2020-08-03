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

#include "_braid.h"

#define TAG_DATA_TO_MASTER (1)
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

double my_Residual(braid_App app, double u, int index) { return 0.0; }

double objective(double alpha, double *us_prev, double *us, int len,
                 braid_App app) {
    double result = 0.0;
    for (int i = 0; i < len; i++) {
        double residual =
            my_Residual(app, alpha * us_prev[i] + (1 - alpha) * us[i], i);
        result += residual * residual;
    }
    return result;
}

typedef double (*Objective)(double, double *, double *, int, braid_App);

double line_search(Objective f, double *us_prev, double *us, int len,
                   braid_App app) {
    double a = 0.0;
    double b = 1.0;
    double gr = (1.0 + sqrt(5.0)) / 2.0;
    double c = 1 - 1 / gr;
    double d = 1 / gr;
    while (fabs(c - d) > 0.00001) {
        if (f(c, us_prev, us, len, app) < f(d, us_prev, us, len, app)) {
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
    int num_vectors;
    int level;
    // The lower and upper indices in time that our proccess has computed
    int lower_t, upper_t;

    braid_SyncStatusGetLevel(status, &level);
    braid_SyncStatusGetNTPoints(status, &num_vectors);
    
    if (level != 0) {
        return 0;
    }

    braid_SyncStatusGetTIUL(status, &upper_t, &lower_t, level);

    braid_Vector *us = malloc(sizeof(braid_Vector) * (upper_t - lower_t));

    _braid_Grid **grids = _braid_StatusElt(status, grids);

    _braid_Grid *grid = grids[level];

    braid_BaseVector *almost = grid->ua;

    for (int i = 0; i < upper_t - lower_t; i++) {
        us[i] = almost[i]->userVector;
    }

    if (app->rank == 0) {
        // Master process
        // Receive all other vectors
        // This code will be run on the first braid_Vector
        // belonging to process 0
    } else {
        // Worker process
        // Send u's data to master process
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

    /* For our custom line searching */
    braid_SetSync(core, my_Sync);

    /* Set some typical Braid parameters */
    braid_SetPrintLevel(core, 2);
    braid_SetMaxLevels(core, 2);
    braid_SetAbsTol(core, 1.0e-06);
    braid_SetCFactor(core, -1, 2);

    /* Run simulation, and then clean up */
    braid_Drive(core);

    braid_Destroy(core);
    free(app);
    MPI_Finalize();

    return (0);
}
