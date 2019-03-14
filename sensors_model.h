#ifndef SENSORS_MODEL_H_
#define SENSORS_MODEL_H_

#include "utils.h"

// Types and structs
typedef struct
{
    double lat_lon_noise_std_dev;
    double alt_noise_std_dev;
    double speed_noise_std_dev;

    double lat_lon_alt[3]; // deg, deg, m
    double gps_speed[3]; // m/s
    double ground_speed; // m/s
    double eph; // GpsSensor HDOP (m)
    double epv; // GpsSensor VDOP (m)
    int fix; // 0-1: No fix ; 2: 2D fix ; 3: 3D fix
    double cog; // course over ground (deg)
    int visible_sats; // number of visible satellites.
} GpsSensor;

typedef struct
{
    double acc_noise_std_dev;
    double gyro_noise_std_dev;

    double acc[3]; // accelerometer (m/(s^2))
    double gyro[3]; // gyroscope (rad/s)
} ImuSensor;

typedef struct
{
    double mag_noise_std_dev;
    double mag_field[3]; // magnetic field (gauss)
} MagSensor;

typedef struct
{
    double baro_noise_std_dev;

    double temperature; // degC
    double pressure; // barometric pressure (Pa)
    double diff_pressure; // Differential pressure
    double pressure_alt; // Altitude from pressure
} BaroSensor;

typedef struct
{
    GpsSensor gps;
    ImuSensor imu;
    MagSensor mag;
    BaroSensor baro;
} Sensors;

// Functions
void init_gps(GpsSensor *gps, double eph, double epv, double fix, double visible_sats, double lat_lon_noise_std_dev, double alt_noise_std_dev, double speed_noise_std_dev);
void update_gps(GpsSensor *gps, double pos_e[3], double vel_e[3]);

void init_imu(ImuSensor *imu, double acc_noise_std_dev, double gyro_noise_std_dev);
void update_imu(ImuSensor *imu, double acc_b[3], double omega_b[3], double dcm_be[3][3]);

void init_mag(MagSensor *mag, double mag_noise_std_dev);
void update_mag(MagSensor *mag, double lat, double lon, double dcm_be[3][3]);

void init_baro(BaroSensor *baro, double baro_noise_std_dev);
void update_baro(BaroSensor *baro, double alt, double vel_b_x);

#endif