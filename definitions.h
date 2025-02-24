#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  EULER
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            2

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   EIGHT_WAVES
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 YES

/* -- user-defined parameters (labels) -- */

#define  GNUPLOT_HEADER                 YES
#define  VTK_TIME_INFO                  YES
#define  VTK_VECTOR_DUMP                YES

#define  DATESHIFT                      0
#define  DAILYBC                        1
#define  USE_POLARITY                   2

/* [Beg] user-defined constants (do not change this line) */

#define UNIT_LENGTH (CONST_au) // AU
#define UNIT_VELOCITY (CONST_au / 86400.0) // au/day
#define UNIT_DENSITY (CONST_mp)   // assume protons

/* [End] user-defined constants (do not change this line) */
