#ifndef PX4_QUAD_SIM_H_
#define PX4_QUAD_SIM_H_

#include "quad_model.h"

// Constants.
#define HIL 0
#define SIL_SPEED_FACTOR 3

#define SENSOR_FREQ 250
#define GPS_FREQ 50

// Global variables.
Quad quad;

// Function prototypes.
int init_sim();
int advance_sim(uint64_t time_usec, uint64_t prev_time_usec, double wind_vel_e[3]);

#endif