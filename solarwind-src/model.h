#ifndef SOLARWIND_MODEL_H_
#define SOLARWIND_MODEL_H_

#include "pluto.h"

#include "solarwind-src/bnd.h"
#include "solarwind-src/utils.h"

double rotate_bc(double k, const double t) {
  // Carrington sidereal rotation rate
  k -= (180 * t / 25.38);
  // Take into account frame rotation and get synodic rotation rate (27.2753)
  k += (t * 180 * g_OmegaZ / (2 * CONST_PI));
  
  while (k < 0) k += 180;
  while (k >= 180) k -= 180;

  return k;
}

double lerp(double a, double b, double s) {
    return a * (1 - s) + b * s;
}

double normalize_time(double t, double left, double right) {
  return (t - left) / (right - left);
}

void normalize(double *D, double *T, double *B1, double *B3, const double x1) {
  const double coef1 = 0.1 / (x1);
  const double coef2 = 0.1 * 0.1 / (x1 * x1);
  const double coef3 = 0.1 * 0.1 * 0.1 / (x1 * x1 * x1);

  *D *= coef2;
  *T *= coef1;
  *B1 *= coef2;
  *B3 *= coef2;
}

void interpolate_vars(const boundary_data* solarwind_data, const double kk, const int jj, 
                      const int frame, double *D, double *V1, double *T, double *B1, double *B3) {
  int kk0 = (int) kk;
  int kk1 = min(179, kk0 + 1);

  double s = kk - kk0;
  int idx0 = get_data_index(frame, kk0, jj);
  int idx1 = get_data_index(frame, kk1, jj);
  *D  = lerp(solarwind_data->D[idx0],  solarwind_data->D[idx1],  s);
  *V1 = lerp(solarwind_data->V1[idx0], solarwind_data->V1[idx1], s);
  *T  = lerp(solarwind_data->T[idx0],  solarwind_data->T[idx1],  s);
  *B1 = lerp(solarwind_data->B1[idx0], solarwind_data->B1[idx1], s);
  *B3 = lerp(solarwind_data->B3[idx0], solarwind_data->B3[idx1], s);
}

void interpolate_ambient(const boundary_data* solarwind_data, const int k, const int j,
                         const double t, const double x1,
                         double* D, double* V1, double* T, double* B1, double* B3) {
  const double bkg_frame_time = convert_to_pluto_time(solarwind_data->TIME[solarwind_data->bkg_frame]);
	const double kk = rotate_bc(k, -bkg_frame_time + t);

  int kk0 = (int) kk;
  int kk1 = min(179, kk0 + 1);

  double s = kk - kk0;
  int idx0 = get_data_index(solarwind_data->bkg_frame, kk0, j);
  int idx1 = get_data_index(solarwind_data->bkg_frame, kk1, j);
  *D  = lerp(solarwind_data->D[idx0],  solarwind_data->D[idx1],  s);
  *V1 = lerp(solarwind_data->V1[idx0], solarwind_data->V1[idx1], s);
  *T  = lerp(solarwind_data->T[idx0],  solarwind_data->T[idx1],  s);
  *B1 = lerp(solarwind_data->B1[idx0], solarwind_data->B1[idx1], s);
  *B3 = lerp(solarwind_data->B3[idx0], solarwind_data->B3[idx1], s);
  normalize(D, T, B1, B3, x1);
}

void interpolate_cme(const boundary_data* solarwind_data, const int left_frame, const int right_frame,
                     const int k, const int j, const double t, const double x1,
                     double* D, double* V1, double* T, double* B1, double* B3) {
  const double left_time = convert_to_pluto_time(solarwind_data->TIME[left_frame]);
  const double right_time = convert_to_pluto_time(solarwind_data->TIME[right_frame]);

  const double kk_cur = rotate_bc(k, t - left_time);
  const double kk_next = rotate_bc(k, -right_time + t);

  double Dcur, V1cur, Tcur, B1cur, B3cur;
  interpolate_vars(solarwind_data, kk_cur, j, left_frame, &Dcur, &V1cur, &Tcur, &B1cur, &B3cur);
  double Dnext, V1next, Tnext, B1next, B3next;
  interpolate_vars(solarwind_data, kk_next, j, right_frame, &Dnext, &V1next, &Tnext, &B1next, &B3next);

  double q = normalize_time(t, left_time, right_time);
  if (fabs(q) > 1) {
      printLog("[INTERP COEF ERROR interpolate_cme] q: %lf, t: %lf, left: %lf, right: %lf",
                q, t, left_time, right_time);
      exit(-1);
  }

  *D  = lerp(Dcur,  Dnext,  q);
  *V1 = lerp(V1cur, V1next, q);
  *T  = lerp(Tcur,  Tnext,  q);
  *B1 = lerp(B1cur, B1next, q);
  *B3 = lerp(B3cur, B3next, q);
  normalize(D, T, B1, B3, x1);
}

#endif