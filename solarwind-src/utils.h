#ifndef SOLARWIND_UTILS_H_
#define SOLARWIND_UTILS_H_

#include "pluto.h"

double convert_to_pluto_time(double time) {
    return time / 86400.0;
}

int get_data_index(int frame, int k, int j) {
  return (frame * 180 + k) * 60 + j;
}

int min(int x, int y) {
    return x < y ? x : y;
}

int max(int x, int y) {
    return x < y ? y : x;
}

int obsdate_cal_to_seconds_from_midnight(const char* date) {
    char hour[3] = {date[11], date[12], '\0'};
    return atoi(hour) * 60 * 60;
}

void map_to_global_indexes(const int local_j, const int local_k, double* global_j, double* global_k) {
    int total_pieces_k = 180 / NX3;
    *global_j = local_j - 2 + prank / total_pieces_k * NX2;
    *global_k = local_k - 2 + prank % total_pieces_k * NX3;
    if (*global_j < 0)      *global_j = 0;
    if (*global_j > 59)     *global_j = 59;
    while (*global_k < 0)   *global_k += 180;
    while (*global_k > 179) *global_k -= 180;
}

#endif
