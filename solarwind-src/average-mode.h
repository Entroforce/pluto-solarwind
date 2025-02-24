#ifndef SOLARWIND_AVERAGE_MODE_H_
#define SOLARWIND_AVERAGE_MODE_H_

#include "pluto.h"

#include "solarwind-src/cme-timeline.h"
#include "solarwind-src/utils.h"

void average_boundary(boundary_data* today_solarwind_data, boundary_data** daily_solarwind_data,
                      cme_timeline_t* cme_timeline,
                      Grid* grid, int local_j, int local_k,
                      double* D, double* V1, double* T, double* B1, double* B3, double x1, int step_count) {
    double j, k;
    map_to_global_indexes(local_j, local_k, &j, &k);

    double t = g_time + g_inputParam[DATESHIFT];
    int cme_index = get_cme_index_by_pluto_time(cme_timeline, t);

    if (cme_index >= 0 && cme_index < (int) cme_timeline->len) {
        const cme_segment_t* cme_segment = &cme_timeline->cme_segments[cme_index];
        const boundary_data* cme_solarwind_data = daily_solarwind_data[cme_segment->daily_idx];

        if (step_count == MAX_DEBUG_STEPS) {
            printLog("[DAILY CME] daily_idx: %d, cme_index: %d, cme_left_time: %lf, cme_right_right: %lf, t: %lf, bkg_frame: %d, bkg_frame_time: %lf, ok: %d\n",
                      cme_segment->daily_idx, cme_index, cme_segment->left_time, cme_segment->right_time, t, today_solarwind_data->bkg_frame, today_solarwind_data->TIME[today_solarwind_data->bkg_frame],
                      (today_solarwind_data == cme_solarwind_data));
        }

        interpolate_cme(cme_solarwind_data,
                        cme_segment->left_frame_index,
                        cme_segment->right_frame_index,
                        k, j,
                        t,
                        x1, D, V1, T, B1, B3);
    } else {
        if (step_count == MAX_DEBUG_STEPS) {
            printLog("[AVERAGE AMBIENT] t: %lf\n", t);
        }

        interpolate_ambient(today_solarwind_data, k, j, t, x1, D, V1, T, B1, B3);
    }
}

#endif
