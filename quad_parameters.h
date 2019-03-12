#ifndef QUAD_PARAMETERS_H_
#define QUAD_PARAMETERS_H_

#include <math.h>
#include "utils.h"

// Environment
#define HOME_LAT 47.397742
#define HOME_LON 8.545594
#define HOME_ALT 488.0
#define HOME_YAW 0.0*M_PI/180.0

#define NOISE 0

// Quad physical properties
#define MASS 4.5
#define I_XX 0.23
#define I_YY 0.235
#define I_ZZ 0.328

#define D_ARM 0.49
#define R_D 0.018
#define C_D 0
#define THRUST_TC 0.07

#define THRUST_MAX_FORCE 4 * 3.2 * GRAVITY
#define THRUST_HOVER_NORM (MASS * GRAVITY) / (THRUST_MAX_FORCE)

// Quad sensor parameters
#define GPS_EPH 0.3
#define GPS_EPV 0.4
#define GPS_FIX 3
#define GPS_NUM_SATS 10
#define GPS_LAT_LON_NOISE 1e-6
#define GPS_ALT_NOISE 0.01
#define GPS_SPEED_NOISE 0.01

#define IMU_ACC_NOISE 0.05
#define IMU_GYRO_NOISE 0.01

#define MAG_DECL 2.13
#define MAG_INCL 63.32
#define MAG_SCALE 2
#define MAG_NOISE 0.005

#define BARO_TEMP 32

#endif