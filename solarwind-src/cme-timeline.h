#ifndef SOLARWIND_CME_TIMELINE_H_
#define SOLARWIND_CME_TIMELINE_H_

#include "solarwind-src/bnd.h"
#include "solarwind-src/utils.h"

#define MAX_CME_SEGMENTS 1000

typedef struct {
    double left_time;
    double right_time;
    int daily_idx;
    int left_frame_index;
    int right_frame_index;
} cme_segment_t;

typedef struct {
    size_t len;
    cme_segment_t cme_segments[MAX_CME_SEGMENTS];
} cme_timeline_t;

int cme_segment_comparator(const void* a, const void* b) {
    cme_segment_t* cme_segment_left = (cme_segment_t*) a;
    cme_segment_t* cme_segment_right = (cme_segment_t*) b;

    return cme_segment_left->left_time > cme_segment_right->left_time
        || cme_segment_left->left_time == cme_segment_right->left_time &&
           cme_segment_left->right_time > cme_segment_right->right_time
        || cme_segment_left->left_time == cme_segment_right->left_time &&
           cme_segment_left->right_time == cme_segment_right->right_time &&
           cme_segment_left->daily_idx > cme_segment_right->daily_idx;
}

// НЕ УЧИТЫВАЕТСЯ BKG_FRAME!!!!!!
cme_timeline_t* create_timeline(boundary_data** daily_solarwind_data, size_t boundaries_amount) {
    cme_timeline_t* cme_timeline = (cme_timeline_t*) malloc(sizeof(cme_timeline_t));

    size_t current = 0;
    //// не удалять, нужно на будущее!
    // for (int daily_idx = 0; daily_idx < boundaries_amount; ++daily_idx) {
    //     boundary_data* current_solarwind_data = daily_solarwind_data[daily_idx];
    //     for (size_t i = 0; i < current_solarwind_data->ntime - 1; ++i) { // ТУТ СТАРЫЙ СПОСОБ ОБРАБОТКИ BKG_FRAME
    //         cme_segment_t* cme_segment = &cme_timeline->cme_segments[current++];
    //         cme_segment->left_time = current_solarwind_data->pluto_time_from_main_bnd + convert_to_pluto_time(current_solarwind_data->TIME[i]);
    //         cme_segment->right_time = current_solarwind_data->pluto_time_from_main_bnd + convert_to_pluto_time(current_solarwind_data->TIME[i + 1]);
    //         cme_segment->daily_idx = daily_idx;
    //         cme_segment->left_frame_index = i;
    //         cme_segment->right_frame_index = i + 1;
    //     }
    // }

    for (int daily_idx = 0; daily_idx < boundaries_amount; ++daily_idx) {
        boundary_data* current_solarwind_data = daily_solarwind_data[daily_idx];
        int left = 0, right = 1;
        while (right < current_solarwind_data->ntime) {
            if (right == current_solarwind_data->bkg_frame) {
                ++right;
            }
            if (right >= current_solarwind_data->ntime) {
                break;
            }

            cme_segment_t* cme_segment = &cme_timeline->cme_segments[current++];
            cme_segment->left_time = current_solarwind_data->pluto_time_from_main_bnd + convert_to_pluto_time(current_solarwind_data->TIME[left]);
            cme_segment->right_time = current_solarwind_data->pluto_time_from_main_bnd + convert_to_pluto_time(current_solarwind_data->TIME[right]);
            cme_segment->daily_idx = daily_idx;
            cme_segment->left_frame_index = left;
            cme_segment->right_frame_index = right;

            left = right; ++right;
        }
    }

    cme_timeline->len = current;

    qsort(cme_timeline->cme_segments, cme_timeline->len, sizeof(cme_segment_t), cme_segment_comparator);
    return cme_timeline;
}

int get_cme_index_by_pluto_time(cme_timeline_t* cme_timeline, const double t) {
    int cme_index = -1;

    for (int i = 0; i < (int) cme_timeline->len; ++i) {
        if (t >= cme_timeline->cme_segments[i].left_time && t <= cme_timeline->cme_segments[i].right_time)
        {
            cme_index = i;
            break;
        }
    }

    return cme_index;
}

#endif