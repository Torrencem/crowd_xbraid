
#include "_braid.h"
#include "util.h"
#include <assert.h>

#define TAG_US_PREV 4

#ifndef LINE_SEARCH_A_START
#define LINE_SEARCH_A_START 0.0
#endif

#ifndef LINE_SEARCH_B_START
#define LINE_SEARCH_B_START 2.0
#endif

#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)

double compute_total_residual(braid_Core core, int level) {
    braid_App app = _braid_CoreElt(core, app);
    _braid_Grid **grids = _braid_CoreElt(core, grids);
    braid_Real *ta = _braid_GridElt(grids[level], ta);
    braid_Int ilower = _braid_GridElt(grids[level], ilower);
    braid_Int iupper = _braid_GridElt(grids[level], iupper);

    double result = 0.0;
    for (int i = ilower; i <= iupper; i++) {
        braid_BaseVector r;
        _braid_BaseInit(core, app, ta[i], &r);
        _braid_TriResidual(core, level, i, 1, &r);
        double rnorm;
        _braid_BaseSpatialNorm(core, app, r, &rnorm);
        result += rnorm;
    }

    return result;
}

int line_search_sync(braid_App app, braid_SyncStatus status) {
    MPI_Status mpi_status;
    MPI_Request request;
    int level;
    // The lower and upper indices in time that our proccess has computed
    int lower_t, upper_t;
    int num_procesors;
    int num_total_points;
    int iter;
    braid_SyncStatusGetIter(status, &iter);
    braid_Core core = (braid_Core)status;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procesors);

    braid_SyncStatusGetLevel(status, &level);

    braid_SyncStatusGetNTPoints(status, &num_total_points);

    if (level != 0) {
        return 0;
    }

    // Collect us - the vectors for our process
    _braid_Grid **grids = _braid_StatusElt(status, grids);
    _braid_Grid *grid = grids[level];
    upper_t = grid->iupper;
    lower_t = grid->ilower;

    int bvector_size;
    _braid_BaseBufSize(core, app, &bvector_size, (braid_BufferStatus)core);

    int my_num_values = upper_t - lower_t + 1; // Assuming upper_t is inclusive

    braid_BaseVector *us = malloc(my_num_values * sizeof(braid_BaseVector));

    for (int i = 0; i < my_num_values; i++) {
        braid_BaseVector base_v;

        _braid_UGetVectorRef((braid_Core)status, level, lower_t + i, &base_v);

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
        // https://en.wikipedia.org/wiki/Golden-section_search
        double tol = 1e-2;
        double a = min(LINE_SEARCH_A_START, LINE_SEARCH_B_START);
        double b = max(LINE_SEARCH_A_START, LINE_SEARCH_B_START);
        double h = b - a;
        assert(h >= tol);
        double invphi = (-1.0 + sqrt(5.0)) / 2.0;
        double invphi2 = (3.0 - sqrt(5.0)) / 2.0;
        int n = ceil(log(tol / h) / log(invphi));

        double c = a + invphi2 * h;
        double d = a + invphi * h;

        // Update our us values to be with alpha = c
        for (int i = 0; i < my_num_values; i++) {
            _braid_BaseSum(core, app, 1.0, us_fut[i], 0.0, us[i]);
            _braid_BaseSum(core, app, 1.0 - c, us_prev[i], c, us[i]);
        }
        // Compute our part of the residual
        double my_res_with_c = compute_total_residual(core, level);
        double yc;
        MPI_Allreduce(&my_res_with_c, &yc, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        
        // Update our us values to be with alpha = d
        for (int i = 0; i < my_num_values; i++) {
            _braid_BaseSum(core, app, 1.0, us_fut[i], 0.0, us[i]);
            _braid_BaseSum(core, app, 1.0 - d, us_prev[i], d, us[i]);
        }
        // Compute our part of the residual
        double my_res_with_d = compute_total_residual(core, level);
        double yd;
        MPI_Allreduce(&my_res_with_d, &yd, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

        for (int k = 0; k < n - 1; k++) {
            if (yc < yd) {
                b = d;
                d = c;
                yd = yc;
                h = invphi * h;
                c = a + invphi2 * h;
                // Update our us values to be with alpha = c
                for (int i = 0; i < my_num_values; i++) {
                    _braid_BaseSum(core, app, 1.0, us_fut[i], 0.0, us[i]);
                    _braid_BaseSum(core, app, 1.0 - c, us_prev[i], c, us[i]);
                }
                // Compute our part of the residual
                double my_res_with_c = compute_total_residual(core, level);
                MPI_Allreduce(&my_res_with_c, &yc, 1, MPI_DOUBLE, MPI_SUM,
                              MPI_COMM_WORLD);
            } else {
                a = c;
                c = d;
                yd = yc;
                h = invphi * h;
                d = a + invphi * h;
                // Update our us values to be with alpha = d
                for (int i = 0; i < my_num_values; i++) {
                    _braid_BaseSum(core, app, 1.0, us_fut[i], 0.0, us[i]);
                    _braid_BaseSum(core, app, 1.0 - d, us_prev[i], d, us[i]);
                }
                // Compute our part of the residual
                double my_res_with_d = compute_total_residual(core, level);
                MPI_Allreduce(&my_res_with_d, &yd, 1, MPI_DOUBLE, MPI_SUM,
                              MPI_COMM_WORLD);
            }
        }
        
        double alpha = yc < yd ? (a + d) / 2 : (c + b) / 2;
        if (rank == 0)
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
    // free(us_prev_message);

    return 0;
}
