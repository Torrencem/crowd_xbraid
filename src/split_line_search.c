
#include "_braid.h"
#include "util.h"
#include <assert.h>

#define TAG_US_PREV 4

double compute_total_residual(braid_Core core, int level) {
    MPI_Comm comm = _braid_CoreElt(core, comm);
    braid_App app = _braid_CoreElt(core, app);
    braid_Real tol = _braid_CoreElt(core, tol);
    braid_Int iter = _braid_CoreElt(core, niter);
    _braid_Grid **grids = _braid_CoreElt(core, grids);
    braid_StepStatus status = (braid_StepStatus)core;
    braid_Int nrefine = _braid_CoreElt(core, nrefine);
    braid_Int ncpoints = _braid_GridElt(grids[level], ncpoints);
    braid_Int gupper = _braid_CoreElt(core, gupper);
    braid_Int tnorm = _braid_CoreElt(core, tnorm);
    braid_Real *ta = _braid_GridElt(grids[level], ta);
    braid_Int ilower = _braid_GridElt(grids[level], ilower);
    _braid_CommHandle *send_handle;
    braid_Int send_index;

    braid_Int flo, fhi, fi, ci, ii, interval;
    braid_Real rnorm_temp, rnorm = 0, global_rnorm = 0;
    braid_BaseVector u, r;

    _braid_UCommInit(core, level);

    /* Start from the right-most interval. */
    for (interval = ncpoints; interval > -1; interval--) {
        _braid_GetInterval(core, level, interval, &flo, &fhi, &ci);
        if (flo <= fhi) {
            _braid_UGetVector(core, level, flo - 1, &u);
        } else if (ci > _braid_CoreElt(core, initiali)) {
            _braid_UGetVector(core, level, ci - 1, &u);
        }

        /* Generate F-points and get residual. */
        for (fi = flo; fi <= fhi; fi++) {
            _braid_BaseClone(core, app, u, &r);
            _braid_Step(core, level, fi, NULL, u);

            /* Update local processor norm. */
            ii = fi - ilower;
            _braid_StepStatusInit(ta[ii - 1], ta[ii], fi - 1, tol, iter, level,
                                  nrefine, gupper, status);
            /* _braid_BaseFullResidual(core, app, r, u, status); */
            /* _braid_Residual(core, level, ii, u, r); */
            // From above
            //

            braid_BaseVector rstop;
            if (_braid_CoreElt(core, residual) == NULL) {
                /* By default: r = ustop - \Phi(ustart)*/
                _braid_GetUInit(core, level, ii, r, &rstop);
                _braid_BaseStep(core, app, rstop, NULL, r, level, status);
                _braid_BaseSum(core, app, 1.0, u, -1.0, r);
            } else {
                /* Call the user's residual routine */
                _braid_BaseResidual(core, app, u, r, status);
            }

            //
            _braid_BaseSpatialNorm(core, app, r, &rnorm_temp);
            if (tnorm == 1) /* one-norm */
            {
                rnorm += rnorm_temp;
            } else if (tnorm == 3) /* inf-norm */
            {
                rnorm = (((rnorm_temp) > (rnorm)) ? (rnorm_temp) : (rnorm));
            } else /* default two-norm */
            {
                rnorm += (rnorm_temp * rnorm_temp);
            }

            /* Communicate w/ neighbor nodes. */
            send_handle = _braid_GridElt(grids[level], send_handle);
            send_index = _braid_GridElt(grids[level], send_index);
            if (fi == send_index) {
                /* Post send to neighbor processor */
                _braid_CommSendInit(core, level, fi, u, &send_handle);
                _braid_GridElt(grids[level], send_index) = _braid_SendIndexNull;
                _braid_GridElt(grids[level], send_handle) = send_handle;
            }
            _braid_BaseFree(core, app, r);
        }
        /* Residual from C-point. */
        if (ci > _braid_CoreElt(core, initiali)) {
            /* Update local processor norm. */
            ii = ci - ilower;
            _braid_StepStatusInit(ta[ii - 1], ta[ii], ci - 1, tol, iter, level,
                                  nrefine, gupper, status);
            _braid_UGetVector(core, level, ci, &r);
            /* _braid_BaseFullResidual(core, app, r, u, status); */
            /* _braid_Residual(core, level, ii, u, r); */
            // From above
            //

            braid_BaseVector rstop;
            if (_braid_CoreElt(core, residual) == NULL) {
                /* By default: r = ustop - \Phi(ustart)*/
                _braid_GetUInit(core, level, ii, r, &rstop);
                _braid_BaseStep(core, app, rstop, NULL, r, level, status);
                _braid_BaseSum(core, app, 1.0, u, -1.0, r);
            } else {
                /* Call the user's residual routine */
                _braid_BaseResidual(core, app, u, r, status);
            }

            //
            _braid_BaseSpatialNorm(core, app, u, &rnorm_temp);

            if (tnorm == 1) /* one-norm */
            {
                rnorm += rnorm_temp;
            } else if (tnorm == 3) /* inf-norm */
            {
                rnorm = (((rnorm_temp) > (rnorm)) ? (rnorm_temp) : (rnorm));
            } else /* default two-norm */
            {
                rnorm += (rnorm_temp * rnorm_temp);
            }
            _braid_BaseFree(core, app, r);
        }

        if ((flo <= fhi) || (ci > _braid_CoreElt(core, initiali))) {
            _braid_BaseFree(core, app, u);
        }
    }

    _braid_UCommWait(core, level);

    /* Compute global residual norm. */
    if (tnorm == 1) /* one-norm reduction */
    {
        MPI_Allreduce(&rnorm, &global_rnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
    } else if (tnorm == 3) /* inf-norm reduction */
    {
        MPI_Allreduce(&rnorm, &global_rnorm, 1, braid_MPI_REAL, MPI_MAX, comm);
    } else /* default two-norm reduction */
    {
        MPI_Allreduce(&rnorm, &global_rnorm, 1, braid_MPI_REAL, MPI_SUM, comm);
        global_rnorm = sqrt(global_rnorm);
    }

    return global_rnorm;
}

int line_search_sync(braid_App app, braid_SyncStatus status) {
    MPI_Status mpi_status;
    MPI_Request request;
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

    // Collect us - the vectors for our process
    _braid_Grid **grids = _braid_StatusElt(status, grids);
    _braid_Grid *grid = grids[level];
    upper_t = grid->cupper / c_factor;
    lower_t = grid->clower / c_factor;

    int bvector_size;
    _braid_BaseBufSize(core, app, &bvector_size, (braid_BufferStatus)core);

    int my_num_values = upper_t - lower_t + 1; // Assuming upper_t is inclusive

    braid_BaseVector *us = malloc(my_num_values * sizeof(braid_BaseVector));

    for (int i = 0; i < my_num_values; i++) {
        int fine_index;
        _braid_MapCoarseToFine(lower_t + i, c_factor, fine_index);

        braid_BaseVector base_v;

        _braid_UGetVectorRef((braid_Core)status, level, fine_index, &base_v);

        assert(base_v != NULL);
        assert(base_v->userVector != NULL);

        /* base_v->userVector->value = 1000.0; */

        us[i] = malloc(sizeof(struct _braid_BaseVector_struct));

        us[i]->userVector = base_v->userVector;
        us[i]->bar = base_v->bar;
    }

    // Now that we have us...
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (iter != 0) {
        // Receive us_prev from last iteration
        char *us_prev_message = malloc(bvector_size * my_num_values);
        MPI_Recv(us_prev_message, bvector_size * my_num_values, MPI_BYTE, rank,
                 TAG_US_PREV, MPI_COMM_WORLD, &mpi_status);

        // Unpack us_prev from us_prev_message
        braid_BaseVector *us_prev =
            malloc(sizeof(braid_BaseVector) * my_num_values);
        for (int i = 0; i < my_num_values; i++) {
            _braid_BaseBufUnpack(core, app, us_prev_message + i * bvector_size,
                                 us_prev + i, (braid_BufferStatus)core);
        }
        free(us_prev_message);

        // Make a copy of our u values into us_fut
        braid_BaseVector *us_fut =
            malloc(sizeof(braid_BaseVector) * my_num_values);
        for (int i = 0; i < my_num_values; i++) {
            _braid_BaseInit(core, app, grids[level + 1]->ta[i], us_fut + i);
            _braid_BaseClone(core, app, us[i], us_fut + i);
        }

        // Perform the line search in every processor
        double a = 0.0;
        double b = 1.0;
        double gr = (1.0 + sqrt(5.0)) / 2.0;
        double c = b - (b - a) / gr;
        double d = a + (b - a) / gr;

        while (fabs(c - d) > 0.00001) {
            // Update our us values to be with alpha = c
            for (int i = 0; i < my_num_values; i++) {
                _braid_BaseSum(core, app, 1.0, us_fut[i], 0.0, us[i]);
                _braid_BaseSum(core, app, 1.0 - c, us_prev[i], c, us[i]);
            }
            // Compute our part of the residual
            double my_res_with_c = compute_total_residual(core, level);

            // Update our us values to be with alpha = d
            for (int i = 0; i < my_num_values; i++) {
                _braid_BaseSum(core, app, 1.0, us_fut[i], 0.0, us[i]);
                _braid_BaseSum(core, app, 1.0 - d, us_prev[i], d, us[i]);
            }
            // Compute our part of the residual
            double my_res_with_d = compute_total_residual(core, level);

            double res_with_c;
            MPI_Allreduce(&my_res_with_c, &res_with_c, 1, MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);
            double res_with_d;
            MPI_Allreduce(&my_res_with_d, &res_with_d, 1, MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);

            printf("Comparing %f with %f\n", res_with_c, res_with_d);
            if (res_with_c < res_with_d) {
                b = d;
            } else {
                a = c;
            }
            c = b - (b - a) / gr;
            d = a + (b - a) / gr;
        }
        double alpha = (a + b) / 2;
        printf("alpha = %f\n", alpha);
        // Update our us values to correspond to the alpha
        for (int i = 0; i < my_num_values; i++) {
            _braid_BaseSum(core, app, 1.0, us_fut[i], 0.0, us[i]);
            _braid_BaseSum(core, app, 1.0 - alpha, us_prev[i], alpha, us[i]);
        }

        // Free all of us_fut
        for (int i = 0; i < my_num_values; i++) {
            _braid_BaseFree(core, app, us_fut[i]);
        }
        free(us_fut);

        // Free all of us_prev
        for (int i = 0; i < my_num_values; i++) {
            _braid_BaseFree(core, app, us_prev[i]);
        }
        free(us_prev);
    }

    // Send our us so they can be used for us_prev of the next iteration
    char *us_prev_message = malloc(bvector_size * my_num_values);
    // Pack us_prev_message from us
    for (int i = 0; i < my_num_values; i++) {
        _braid_BaseBufPack(core, app, us[i], us_prev_message + i * bvector_size,
                           (braid_BufferStatus)core);
        free(us[i]);
    }
    MPI_Isend(us_prev_message, bvector_size * my_num_values, MPI_BYTE, rank,
              TAG_US_PREV, MPI_COMM_WORLD, &request);

    return 0;
}
