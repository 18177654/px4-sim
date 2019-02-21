#ifndef UTILS_H_
#define UTILS_H_

#include <inttypes.h>

#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define LOG_MSG(M, ...)  printf((ANSI_COLOR_MAGENTA M ANSI_COLOR_RESET), ##__VA_ARGS__)

// Defines
typedef double (*functiontype)(double[], int);

// Constants
#define GRAVITY 9.80665 // gravity
#define Pb 101325.0 // static pressure at sea level [Pa]
#define Tb 288.15 // standard temperature at sea level [K]
#define Lb -0.0065 // standard temperature lapse rate [K/m]
#define M 0.0289644 // molar mass of Earth's air [kg/mol]
#define R 8.31432 // universal gas constant
#define ETA_W 0.2 // eta
#define RHO 1.225 // air density
#define R_EARTH 6371000.0 // radius of earth (6378137.0)

// Utility functions
// Prototypes
double deg2rad(double deg);
double rad2deg(double rad);
double zero_mean_noise(double std_dev);
double rand_gauss(double mean, double std_dev);

void ned_to_latlonalt(double ned[3], double latlonalt[3], double home_lat, double home_lon, double home_alt);

void quat_to_euler(double q[4], double euler[3]);
void euler_to_quat(double euler[3], double quat[4]);

void calc_dcm_be(double q[4], double dcm_be[3][3]);
void body_to_earth_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3]);
void earth_to_body_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3]);

double wrap_angle_2pi(double angle);
double wrap_angle_pi(double angle);

void integrate_rk4(functiontype derivs[], double states[], double integ[], int len, double h);
double integrate_euler(double sum, double val, double dt);

uint64_t get_time_usec();

#endif