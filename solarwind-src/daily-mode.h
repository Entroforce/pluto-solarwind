#ifndef SOLARWIND_DAILY_MODE_H_
#define SOLARWIND_DAILY_MODE_H_

#include "pluto.h"

#include "solarwind-src/bnd.h"
#include "solarwind-src/utils.h"
#include "solarwind-src/cme-timeline.h"

static void shift_time_relative_to_main_bnd(
        boundary_data* daily_solarwind, boundary_data* today_solarwind,
        int days_diff) {
    int obsdate_hour_diff = today_solarwind->obsdate_hour - daily_solarwind->obsdate_hour;
    for (size_t i = 0; i < daily_solarwind->ntime; ++i) {
        daily_solarwind->pluto_time_from_main_bnd = convert_to_pluto_time(-(obsdate_hour_diff * 60 * 60 + days_diff * 24 * 60 * 60));
    }
}

boundary_data** read_daily_data(boundary_data* today_solarwind_data) {
    const int daily_bnds_count = -g_inputParam[DATESHIFT];
    char filepath[200] = {0};

    boundary_data** daily_solarwind_data = (boundary_data**) malloc((daily_bnds_count + 1) * sizeof(boundary_data*));
    for (int daily_idx = 1; daily_idx <= daily_bnds_count; ++daily_idx) {
        snprintf(filepath, sizeof(filepath), "./bnds/bnd-%d.nc", daily_idx);
        daily_solarwind_data[daily_idx] = read_bnd(filepath);
        shift_time_relative_to_main_bnd(daily_solarwind_data[daily_idx], today_solarwind_data, daily_idx);
    }

    return daily_solarwind_data;
}


int get_daily_idx_by_time(boundary_data** daily_solarwind_data, const double t) {
    for (int i = 1; i <= 10; ++i) {
        if (t > daily_solarwind_data[i]->pluto_time_from_main_bnd) {
            return i;
        }
    }

    return -1;
}

void process_daily_cme(
        boundary_data** daily_solarwind_data, cme_timeline_t* cme_timeline,
        const int cme_index, const double t, const double daily_solarwind_time,
        const double j, const double k,
        double* D, double* V1, double* T, double* B1, double* B3, double x1) {
    const cme_segment_t* cme_segment = &cme_timeline->cme_segments[cme_index];
    const boundary_data* cme_solarwind_data = daily_solarwind_data[cme_segment->daily_idx];
    const double t_between_cmes
        = cme_solarwind_data->TIME[cme_segment->left_frame_index] + (t - cme_segment->left_time);

    interpolate_cme(cme_solarwind_data,
                    cme_segment->left_frame_index,
                    cme_segment->right_frame_index,
                    k, j,
                    convert_to_pluto_time(t_between_cmes),
                    x1, D, V1, T, B1, B3);
}

void process_daily_ambient(
        boundary_data** daily_solarwind_data, const int daily_idx,
        const double daily_solarwind_time, const double j, const double k,
        double* D, double* V1, double* T, double* B1, double* B3, double x1,
        int force_single_ambient) {
    interpolate_ambient(daily_solarwind_data[daily_idx], k, j, daily_solarwind_time, x1, D, V1, T, B1, B3);

    if (daily_idx > 0 && !force_single_ambient) {
        double Dnext, V1next, Tnext, B1next, B3next;

        interpolate_ambient(daily_solarwind_data[daily_idx - 1], k, j, daily_solarwind_time - 1, x1, &Dnext, &V1next, &Tnext, &B1next, &B3next);

        double q = normalize_time(daily_solarwind_time, 0.0, daily_solarwind_data[daily_idx - 1]->pluto_time_from_main_bnd - daily_solarwind_data[daily_idx]->pluto_time_from_main_bnd);
        if (fabs(q) > 1) {
            printLog("[INTERP COEF ERROR] q: %lf, daily_solarwind_time: %lf, left: 0.0, right: %lf",
                        q, daily_solarwind_time, daily_solarwind_data[daily_idx - 1]->pluto_time_from_main_bnd - daily_solarwind_data[daily_idx]->pluto_time_from_main_bnd);
            exit(-1);
        }                                                                            
        *D  = lerp(*D,  Dnext,  q);
        *V1 = lerp(*V1, V1next, q);
        *T  = lerp(*T,  Tnext,  q);
        *B1 = lerp(*B1, B1next, q);
        *B3 = lerp(*B3, B3next, q);
    }
}

//// !!! может произойти такое, что время самого раннего bnd будет >-10, это надо дополнительно учесть
void daily_boundary(boundary_data** daily_solarwind_data, cme_timeline_t* cme_timeline,
                    Grid* grid, int _i, int local_j, int local_k,
                    double* D, double* V1, double* T, double* B1, double* B3, double x1, int step_count) {
    double j, k;
    map_to_global_indexes(local_j, local_k, &j, &k);

    double t = g_time + g_inputParam[DATESHIFT];
    int daily_idx = get_daily_idx_by_time(daily_solarwind_data, t);

    double daily_solarwind_time;
    int force_single_ambient = 0;
    // если daily_idx == -1, то obsdate самого раннего bnd раньше bnd.nc больше, чем на -10.000
    if (daily_idx == -1) {
        daily_solarwind_time = t - daily_solarwind_data[10]->pluto_time_from_main_bnd;
        force_single_ambient = 1;
        daily_idx = 10; // очень опасно
    } else {
        daily_solarwind_time = t - daily_solarwind_data[daily_idx]->pluto_time_from_main_bnd; // >0!
    }

    int cme_index = get_cme_index_by_pluto_time(cme_timeline, t);
    if (cme_index >= 0 && cme_index < (int) cme_timeline->len) {
        process_daily_cme(daily_solarwind_data, cme_timeline,
                          cme_index, t, daily_solarwind_time,
                          j, k, D, V1, T, B1, B3, x1);
    } else {
        process_daily_ambient(daily_solarwind_data, daily_idx, daily_solarwind_time,
                              j, k, D, V1, T, B1, B3, x1, force_single_ambient);
    }
}

#endif
