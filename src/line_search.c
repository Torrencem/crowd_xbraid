
#include "_braid.h"
#include <assert.h>

#define TAG_DATA_TO_MASTER 2
#define TAG_LINE_SEARCH_RESULT 3
#define TAG_US_PREV 4

double objective(double alpha, braid_BaseVector *us_prev, braid_BaseVector *us,
                 int len, int level, braid_App app, braid_Core core) {
    printf("us (len %d) is: ", len);
    for (int i = 0; i < len; i++) {
        printf("%f, ", us[i]->userVector->value);
    }
    printf("\nus_prev is: ");
    for (int i = 0; i < len; i++) {
        printf("%f, ", us_prev[i]->userVector->value);
    }
    printf("\n");
    int cfactor;
    _braid_GetCFactor(core, level, &cfactor);
    double result = 0.0;
    // r is reused on each iteration. It stores the residual calculation
    struct _braid_BaseVector_struct *r;
    _braid_BaseInit(core, app, 0.0, &r);
    for (int i = 1; i < len; i++) {
        struct _braid_BaseVector_struct *us_updated;
        _braid_BaseClone(core, app, us[i], &us_updated);
        _braid_BaseSum(core, app, 1 - alpha, us_prev[i], alpha, us_updated);
        struct _braid_BaseVector_struct *us_im1_updated;
        _braid_BaseClone(core, app, us[i - 1], &us_im1_updated);
        _braid_BaseSum(core, app, 1 - alpha, us_prev[i - 1], alpha,
                       us_im1_updated);
        // Level + 1 since we're doing the line search on the coarse grid
        /* _braid_FASResidual(core, level + 1, i, us_updated, r); */
        braid_Real tol = _braid_CoreElt(core, tol);
        braid_Int iter = _braid_CoreElt(core, niter);
        _braid_Grid **grids = _braid_CoreElt(core, grids);
        braid_StepStatus status = (braid_StepStatus)core;
        braid_Int nrefine = _braid_CoreElt(core, nrefine);
        braid_Int gupper = _braid_CoreElt(core, gupper);
        braid_Int ilower = _braid_GridElt(grids[level + 1], ilower);
        braid_Real *ta = _braid_GridElt(grids[level + 1], ta);

        braid_Int ii;

        ii = i - ilower;
        _braid_StepStatusInit(ta[ii - 1], ta[ii], i, tol, iter, level + 2,
                              nrefine, gupper, status);

        // Now compute Phi(us_prev[i])

        braid_BaseVector phi_u_prev;
        _braid_BaseClone(core, app, us_im1_updated, &phi_u_prev);

        /* struct _braid_BaseVector_struct *ustop; */
        /* _braid_BaseClone(core, app, us[i], &ustop); */
        _braid_BaseStep(core, app, NULL, NULL, phi_u_prev, level + 1, status);
        double tmp = phi_u_prev->userVector->value;
        // r = ustop - \Phi(ustart)
        _braid_BaseSum(core, app, -1.0, us_updated, 1.0, phi_u_prev);
        r = phi_u_prev;
        double norm;
        _braid_BaseSpatialNorm(core, app, r, &norm);
        result += norm;
        printf("Residual is %f for index %d (value is %f and %f)\n", norm, i,
               tmp, us_updated->userVector->value);
        /* _braid_BaseFree(core, app, us_updated); */
    }
    printf("Result is %f\n", result);
    return result;
}

typedef double (*Objective)(double, braid_BaseVector *, braid_BaseVector *, int,
                            int, braid_App, braid_Core);

double line_search(Objective f, braid_BaseVector *us_prev, braid_BaseVector *us,
                   int len, int level, braid_App app, braid_Core status) {
    /* return 0.0; */
    double a = 0.0;
    double b = 3.0;
    double gr = (1.0 + sqrt(5.0)) / 2.0;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;
    while (fabs(c - d) > 0.00001) {
        if (f(c, us_prev, us, len, level, app, status) <
            f(d, us_prev, us, len, level, app, status)) {
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

int line_search_sync(braid_App app, braid_SyncStatus status) {
    MPI_Status mpi_status;
    MPI_Request request;
    int num_vectors;
    int level;
    // The lower and upper indices in time that our proccess has computed
    int lower_t, upper_t;
    int num_procesors;
    int num_total_points, c_factor;
    int iter;
    braid_SyncStatusGetIter(status, &iter);
    braid_Core core = (braid_Core)status;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procesors);

    braid_SyncStatusGetLevel(status, &level);

    braid_SyncStatusGetNTPoints(status, &num_total_points);
    _braid_GetCFactor(core, level, &c_factor);

    if (level != 0) {
        return 0;
    }

    num_vectors = num_total_points / c_factor + 1;

    // Collect us - the vectors for our process
    _braid_Grid **grids = _braid_StatusElt(status, grids);
    _braid_Grid *grid = grids[level];
    upper_t = grid->cupper / c_factor;
    lower_t = grid->clower / c_factor;

    int bvector_size;
    _braid_BaseBufSize(core, app, &bvector_size, (braid_BufferStatus)core);
    int safe_message_size =
        (num_vectors / num_procesors + 2) * bvector_size + 2 * sizeof(int);

    int my_num_values = upper_t - lower_t + 1; // Assuming upper_t is inclusive

    assert(2 * sizeof(int) + my_num_values * bvector_size <=
           (unsigned long)safe_message_size);

    braid_BaseVector *us = malloc(my_num_values * sizeof(braid_BaseVector));

    for (int i = 0; i < my_num_values; i++) {
        int fine_index;
        _braid_MapCoarseToFine(lower_t + i, c_factor, fine_index);
        braid_BaseVector base_v;

        _braid_UGetVectorRef((braid_Core)status, level, fine_index, &base_v);

        assert(base_v != NULL);

        us[i] = base_v;
    }

    // Now that we have us...
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        /* printf("I'm the master and I go from %d to %d\n", lower_t, upper_t);
         */
        // Master process
        // Receive all other vectors
        braid_BaseVector *us_combined =
            malloc(sizeof(braid_BaseVector) * num_vectors);
        // handle our own us
        for (int i = lower_t; i <= upper_t; i++) {
            us_combined[lower_t + i] = us[i];
            /* printf("index: %d\n", i); */
        }
        for (int i = 0; i < num_procesors - 1; i++) {
            // Received message:
            // (int: the start index (lower_t)) : (int: the number of us
            // (upper_t - lower_t)) : values (us)
            char *message = malloc(safe_message_size);
            MPI_Recv(message, safe_message_size, MPI_BYTE, MPI_ANY_SOURCE,
                     TAG_DATA_TO_MASTER, MPI_COMM_WORLD, &mpi_status);
            // Finagle the message
            int start_index = ((int *)message)[0];
            int us_len = ((int *)message)[1];
            char *dmessage = (char *)(message + 2 * sizeof(int));
            for (int i = 0; i < us_len; i++) {
                /* us_combined[i + start_index] = dmessage[i]; */
                // use BufUnpack
                int index = i + start_index;
                _braid_BaseBufUnpack(core, app, dmessage + i * bvector_size,
                                     us_combined + index,
                                     (braid_BufferStatus)core);
                /* printf("index: %d\n", index); */
            }
            /* free(dmessage); */
        }
        // Perform the line search
        if (iter != 0) {
            char *us_prev_message = malloc(bvector_size * num_vectors);
            // Receive us_prev from previous iteration of access
            MPI_Recv(us_prev_message, num_vectors * bvector_size, MPI_BYTE, 0,
                     TAG_US_PREV, MPI_COMM_WORLD, &mpi_status);
            // Unpack from message into us_prev
            braid_BaseVector *us_prev =
                malloc(sizeof(braid_BaseVector) * num_vectors);
            for (int i = 0; i < num_vectors; i++) {
                _braid_BaseBufUnpack(core, app,
                                     us_prev_message + i * bvector_size,
                                     us_prev + i, (braid_BufferStatus)core);
            }

            /* free(us_prev_message); */

            double alpha =
                line_search(objective, us_prev, us_combined, num_vectors, level,
                            app, (braid_Core)status);

            /* alpha = 1.0; // TODO Tmp */

            printf("alpha: %f\n", alpha);
            // ... put result in us_combined
            for (int i = 0; i < num_vectors; i++) {
                _braid_BaseSum(core, app, 1 - alpha, us_prev[i], alpha,
                               us_combined[i]);
            }

            // As the main process, we need to update our us[i]'s
            for (int i = lower_t; i <= upper_t; i++) {
                _braid_BaseClone(core, app, us_combined[i], us + i - lower_t);
            }
        }
        // Send us_prev for next iteration of loop
        // Pack us_combined to send
        char *us_combined_packed = malloc(num_vectors * bvector_size);
        for (int i = 0; i < num_vectors; i++) {
            _braid_BaseBufPack(core, app, us_combined[i],
                               us_combined_packed + i * bvector_size,
                               (braid_BufferStatus)core);
        }
        MPI_Isend(us_combined_packed, num_vectors * bvector_size, MPI_BYTE, 0,
                  TAG_US_PREV, MPI_COMM_WORLD, &request);
        // Send updated us_combined to workers
        MPI_Bcast(us_combined_packed, num_vectors * bvector_size, MPI_BYTE, 0,
                  MPI_COMM_WORLD);
        /* free(us_combined); */
        /* free(us_combined_packed); */
    } else {
        /* printf("I'm a worker and I go from %d to %d\n", lower_t, upper_t); */
        // Worker process
        // Send u's data to master process
        char *message = malloc(safe_message_size);
        ((int *)message)[0] = lower_t;
        ((int *)message)[1] = my_num_values;
        char *dmessage = (char *)(message + 2 * sizeof(int));
        for (int i = 0; i < my_num_values; i++) {
            /* dmessage[i] = us[i]->value; */
            // Use BufPack
            _braid_BaseBufPack(core, app, us[i], dmessage + i * bvector_size,
                               (braid_BufferStatus)core);
        }
        MPI_Send(message, safe_message_size, MPI_BYTE, 0, TAG_DATA_TO_MASTER,
                 MPI_COMM_WORLD);
        /* free(message); */
        // Receive updated us values
        char *updated_us = malloc(num_vectors * bvector_size);
        MPI_Bcast(updated_us, num_vectors * bvector_size, MPI_BYTE, 0,
                  MPI_COMM_WORLD);
        // Put them back into us
        for (int i = 0; i < my_num_values; i++) {
            /* us[i]->value = updated_us[i + lower_t]; */
            int index = i + lower_t;
            /* if (index == 20 && iter == 1) { */
            /*     printf("My value is: "); */
            /*     printf("%f", us[i]->userVector->value); */
            /*     printf("\n"); */
            /* } */
            // Use BufUnpack
            _braid_BaseBufUnpack(core, app, updated_us + index * bvector_size,
                                 us + i, (braid_BufferStatus)core);
            /* if (index == 20 && iter == 1) { */
            /*     printf("I got back: "); */
            /*     printf("%f", us[i]->userVector->value); */
            /*     printf("\n"); */
            /* } */
        }
        /* free(updated_us); */
    }

    /* free(us); */

    return 0;
}
