#ifndef SOLARWIND_BND_H_
#define SOLARWIND_BND_H_

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

#include <solarwind-src/utils.h>

#define LATITUDE_DIM 60
#define LONGITUDE_DIM 180

enum BND_MODE {
    AMBIENT,
    CME,
};

typedef struct {
  int left;
  int right;
} segment;

typedef struct {
  int size;
  segment segments[10];
} segments;

typedef struct {
    double *D, *V1, *T, *B1, *B3, *BP, *TIME;
    segments cme_segments;

    int bkg_frame;
    double mean_D, mean_V1, mean_T, mean_B1, mean_B3;

    size_t n2, n3, ntime;
    int mode;

    int obsdate_hour;
    double pluto_time_from_main_bnd;
} boundary_data;

#define CIRCLE_RADIUS 2
int has_circle(double* T, int frame) {
    for (int n2 = CIRCLE_RADIUS; n2 < LATITUDE_DIM - CIRCLE_RADIUS; ++n2) {
        for (int n3 = CIRCLE_RADIUS; n3 < LONGITUDE_DIM - CIRCLE_RADIUS; ++n3) {
            int idx = get_data_index(frame, n3, n2);
            if (T[idx] != 500000.0) { // epsilon
                continue;
            }

            int all_same_horizontal = 0, all_same_vertical = 0;
            for (int radius = 0; radius <= CIRCLE_RADIUS; ++radius) {
                int idx_horizontal = get_data_index(frame, n3 + radius, n2);
                all_same_horizontal = (T[idx_horizontal] == 500000.0); // epsilon
                int idx_vertical = get_data_index(frame, n3, n2 + radius);
                all_same_vertical = (T[idx_vertical] == 500000.0); // epsilon
            }

            if (all_same_horizontal || all_same_vertical) {
                return 1;
            }
        }
    }

    return 0;
}

segments get_cmes_segments(double* T, size_t ntime) {
  segments cmes_segments;

  int* is_cme_frame = malloc(ntime * sizeof(int));
  for (size_t i = 0; i < ntime; ++i) {
    is_cme_frame[i] = has_circle(T, (int) i);
    // printf("is_cme_frame[%d] = %d\n", i, is_cme_frame[i]);
  }

  int left = 0, right = 0, current_segment = 0;
  while (right < ntime) {
    if (is_cme_frame[left] == 0) {
      ++left; ++right;
    } else {
      if (is_cme_frame[right] == 1) {
        ++right;
      } else if (is_cme_frame[right] == 0) {
        segment cme_segment; cme_segment.left = left; cme_segment.right = right - 1;
        --cme_segment.left; ++cme_segment.right; // нужно для учета крайних фреймов, у которых нет "кружка"
        cmes_segments.segments[current_segment++] = cme_segment;
        left = right;
      }
    }
  }

  cmes_segments.size = current_segment;
  return cmes_segments;
}

int get_bkg_frame(const boundary_data* data) {
    int frame = 0;
    for (int i = 0; i < data->cme_segments.size; ++i) {
        if (frame < data->cme_segments.segments[i].left) {
            return frame;
        } else {
            frame = data->cme_segments.segments[i].right + 1;
        }
    }

    if (frame < data->ntime) {
        return frame;
    }

    // mode == CME so cme_segments always has size > 0
    return data->cme_segments.segments[0].right;
}

boundary_data* read_bnd(const char* path) {
    int rc, nc_id;

    rc = nc_open(path, NC_NOWRITE, &nc_id);
    if (rc != NC_NOERR) {
        printf("error opening BC file %s: %d\n", path, rc);
        exit(2);
    }

    boundary_data* data = (boundary_data*) malloc(sizeof(boundary_data));

    int dimension_id, id;
    nc_inq_dimid(nc_id, "n2", &dimension_id);
    nc_inq_dimlen(nc_id, dimension_id, &data->n2);
    nc_inq_dimid(nc_id, "n3", &dimension_id);
    nc_inq_dimlen(nc_id, dimension_id, &data->n3);
    nc_inq_dimid(nc_id, "ntime", &dimension_id);
    nc_inq_dimlen(nc_id, dimension_id, &data->ntime);

    data->D = (double*) malloc(data->ntime * data->n2 * data->n3 * sizeof(double));
    data->V1 = (double*) malloc(data->ntime * data->n2 * data->n3 * sizeof(double));
    data->T = (double*) malloc(data->ntime * data->n2 * data->n3 * sizeof(double));
    data->B1 = (double*) malloc(data->ntime * data->n2 * data->n3 * sizeof(double));
    data->B3 = (double*) malloc(data->ntime * data->n2 * data->n3 * sizeof(double));
    data->BP = (double*) malloc(data->ntime * data->n2 * data->n3 * sizeof(double));
    data->TIME = (double*) malloc(data->ntime * sizeof(double));

    nc_inq_varid(nc_id, "D", &id);
    nc_get_var_double(nc_id, id, data->D);
    nc_inq_varid(nc_id, "V1", &id);
    nc_get_var_double(nc_id, id, data->V1);
    nc_inq_varid(nc_id, "T", &id);
    nc_get_var_double(nc_id, id, data->T);
    nc_inq_varid(nc_id, "B1", &id);
    nc_get_var_double(nc_id, id, data->B1);
    nc_inq_varid(nc_id, "B3", &id);
    nc_get_var_double(nc_id, id, data->B3);
    nc_inq_varid(nc_id, "BP", &id);
    nc_get_var_double(nc_id, id, data->BP);
    nc_inq_varid(nc_id, "TIME", &id);
    nc_get_var_double(nc_id, id, data->TIME);

    if (data->ntime > 1) {
        data->mode = CME;
    } else {
        data->mode = AMBIENT;
    }

	data->bkg_frame = 0;
    if (data->mode == CME) {
        data->cme_segments = get_cmes_segments(data->T, data->ntime);
        data->bkg_frame = get_bkg_frame(data);
    }

    data->mean_D = 0;
    data->mean_V1 = 0;
    data->mean_T = 0;
    data->mean_B1 = 0; 
    data->mean_B3 = 0;

    for (int j = 0; j < data->n2; ++j) {
        for (int k = 0; k < data->n3; ++k) {
            int idx = get_data_index(data->bkg_frame, k, j);

            if (g_inputParam[USE_POLARITY]) {
                data->B1[idx] *= (data->BP[idx] > 0) ? 1 : -1;
                data->B3[idx] *= (data->BP[idx] > 0) ? 1 : -1;
            }

            data->mean_D  += data->D[idx];
            data->mean_V1 += data->V1[idx];
            data->mean_T  += data->T[idx];
            data->mean_B1 += data->B1[idx];
            data->mean_B3 += data->B3[idx];
        }
    }

    data->mean_D /= (data->n2 * data->n3);
    data->mean_V1 /= (data->n2 * data->n3);
    data->mean_T /= (data->n2 * data->n3);
    data->mean_B1 /= (data->n2 * data->n3);
    data->mean_B3 /= (data->n2 * data->n3);

    char obsdate_cal[100];
	nc_get_att_text(nc_id, NC_GLOBAL, "obsdate_cal", obsdate_cal);
    char hour[3] = {obsdate_cal[11], obsdate_cal[12], '\0'};
    data->obsdate_hour = atoi(hour);

    data->pluto_time_from_main_bnd = 0;

    return data;
}

#endif
