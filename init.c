/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include <netcdf.h>
#include "pluto.h"

#include "solarwind-src/bnd.h"
#include "solarwind-src/utils.h"
#include "solarwind-src/model.h"
#include "solarwind-src/daily-mode.h"
#include "solarwind-src/average-mode.h"
#include "solarwind-src/cme-timeline.h"

int step_count = 0;
int bnd_read = 0;
boundary_data* today_solarwind_data;
boundary_data** daily_solarwind_data;
cme_timeline_t* cme_timeline;

void read_bnds() {
	if (bnd_read) {
        return;
    }

    today_solarwind_data = read_bnd("./bnds/bnd.nc");
    if (g_inputParam[DAILYBC]) {
        daily_solarwind_data = read_daily_data(today_solarwind_data);
	    daily_solarwind_data[0] = today_solarwind_data;
        cme_timeline = create_timeline(daily_solarwind_data, 11);
    } else {
        daily_solarwind_data = (boundary_data**) malloc(sizeof(boundary_data*));
	    daily_solarwind_data[0] = today_solarwind_data;
        cme_timeline = create_timeline(&today_solarwind_data, 1);
    }

    for (int i = 0; i < cme_timeline->len; ++i) {
        cme_segment_t* cme_segment = &cme_timeline->cme_segments[i];
        printLog("(%lf, %lf, %d), ", cme_segment->left_time, cme_segment->right_time, cme_segment->daily_idx);
        if (i % 4 == 0) {
            printLog("\n");
        }
    }
    printLog("\n%ld\n", cme_timeline->len);

    bnd_read = 1;
}

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rd dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
    const double coef1 = 0.1 / (x1);
    const double coef2 = 0.1 * 0.1 / (x1 * x1);
    const double coef3 = 0.1 * 0.1 * 0.1 / (x1 * x1 * x1);

    read_bnds();

    v[RHO] = today_solarwind_data->mean_D / CONST_mp / 1000.0 * coef2;
    v[VX1] = today_solarwind_data->mean_V1 / 1000.0 / 149597870.7 * 86400;
    v[VX2] = 0;
    v[VX3] = 0;

    double mu = MeanMolecularWeight(v);
    #if HAVE_ENERGY
    v[PRS] = v[RHO] * (today_solarwind_data->mean_T * coef1) / (KELVIN * mu);
    #endif

    #if PHYSICS == MHD || PHYSICS == RMHD
    v[BX1] = today_solarwind_data->mean_B1 * coef2;
    v[BX2] = 0.0;
    v[BX3] = today_solarwind_data->mean_B3 * coef1;

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = 0.0;
    #endif

    g_OmegaZ = +2.0 * CONST_PI / 365.25;
    g_gamma = 1.5;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
}

#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
    B0[0] = 0.0;
    B0[1] = 0.0;
    B0[2] = 0.0;
}
#endif

/* ********************************************************************ww* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
    const double* x1 = grid->x[IDIR];
    const double* x2 = grid->x[JDIR];
    const double* x3 = grid->x[KDIR];

    int i, j, k;
    if (side == X1_BEG && box->vpos == CENTER) {
        BOX_LOOP(box, k, j, i) {
            double D, V1, T, B1, B3;

            if (g_inputParam[DAILYBC]) {
                daily_boundary(daily_solarwind_data, cme_timeline, grid, i, j, k, &D, &V1, &T, &B1, &B3, x1[i], step_count);
            } else {
                average_boundary(today_solarwind_data, daily_solarwind_data, cme_timeline, grid, j, k, &D, &V1, &T, &B1, &B3, x1[i], step_count);
            }

            if (step_count == MAX_DEBUG_STEPS) {
                step_count = 0;
            }

            d->Vc[RHO][k][j][i] = D / CONST_mp / 1000.0;
            double v = V1 / 1000.0;  // km/s
            d->Vc[VX1][k][j][i] = v / 149597870.7 * 86400; // au/day
            d->Vc[VX2][k][j][i] = 0.;
            // take into account frame rotation
            // (the velocity is radial in the inertial frame but not in rotating frame)
            d->Vc[VX3][k][j][i] = -g_OmegaZ * x1[i] * sin(x2[j]);
            #if HAVE_ENERGY
            double mu = MeanMolecularWeight(NULL);
            d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * T /(KELVIN * mu);
            #endif

            // As per documentation: units must be [Gauss / sqrt(4 pi density velocity^2)]
            // In the bnd.nc, units are teslas, 1 T = 1e4 G
            // see also init.c in Whistler_Waves
            d->Vc[BX1][k][j][i] = B1 * 1e4 / sqrt(4.0 * CONST_PI * UNIT_DENSITY) / UNIT_VELOCITY;
            d->Vc[BX2][k][j][i] = 0.;
            d->Vc[BX3][k][j][i] = B3 * 1e4 / sqrt(4.0 * CONST_PI * UNIT_DENSITY) / UNIT_VELOCITY;
        }
    }

    ++step_count;
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
    double g0 = (UNIT_LENGTH / (UNIT_VELOCITY * UNIT_VELOCITY));
    double r2 = x1 * x1 * UNIT_LENGTH * UNIT_LENGTH;
    g[IDIR] = -(CONST_G * CONST_Msun / r2) * g0;
}

#endif
