#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "sensors_model.h"
#include "quad_parameters.h"

void init_gps(GpsSensor *gps, double eph, double epv, double fix, double visible_sats, double lat_lon_noise_std_dev, double alt_noise_std_dev, double speed_noise_std_dev)
{
    gps->eph = eph;
    gps->epv = epv;
    gps->fix = fix;
    gps->visible_sats = visible_sats;

    gps->lat_lon_noise_std_dev = lat_lon_noise_std_dev;
    gps->alt_noise_std_dev = alt_noise_std_dev;
    gps->speed_noise_std_dev = speed_noise_std_dev;

    gps->lat_lon_alt[0] = HOME_LAT;
    gps->lat_lon_alt[1] = HOME_LON;
    gps->lat_lon_alt[2] = HOME_ALT;
}

void update_gps(GpsSensor *gps, double pos_e[3], double vel_e[3])
{
    int i;

    gps->lat_lon_alt[0] =  (rad2deg((pos_e[0]/R_EARTH)) + HOME_LAT);// + zero_mean_noise(gps->lat_lon_noise_std_dev);
    gps->lat_lon_alt[1] = (rad2deg(pos_e[1]/(R_EARTH*cos(deg2rad(gps->lat_lon_alt[0])))) + HOME_LON);// + zero_mean_noise(gps->lat_lon_noise_std_dev);
    gps->lat_lon_alt[2] = (-pos_e[2] + HOME_ALT);// + zero_mean_noise(gps->alt_noise_std_dev);

    for(i = 0 ; i < 3 ; i++)
        gps->gps_speed[i] = vel_e[i];// + zero_mean_noise(gps->speed_noise_std_dev);

    gps->ground_speed = sqrt(pow(gps->gps_speed[0], 2) + pow(gps->gps_speed[1], 2));
    gps->cog = wrap_angle_2pi(rad2deg(atan2(gps->gps_speed[1], gps->gps_speed[0])));
}

void init_imu(ImuSensor *imu, double acc_noise_std_dev, double gyro_noise_std_dev)
{
    int i;

    imu->acc_noise_std_dev = acc_noise_std_dev;
    imu->gyro_noise_std_dev = gyro_noise_std_dev;

    for(i = 0 ; i < 3 ; i++)
    {
        imu->acc[i] = 0.0f;
        imu->gyro[i] = 0.0f;
    }
}

void update_imu(ImuSensor *imu, double acc_b[3], double omega_b[3], double dcm_be[3][3])
{
    int i;

    for(i = 0 ; i < 3 ; i++)
    {
        imu->acc[i] = acc_b[i] + (dcm_be[i][2] * (-GRAVITY));// + zero_mean_noise(imu->acc_noise_std_dev);
        imu->gyro[i] = omega_b[i];// + zero_mean_noise(imu->gyro_noise_std_dev);
    }
}

void init_mag(MagSensor *mag, double mag_decl, double mag_incl, double mag_scale, double mag_noise_std_dev)
{
    double c_decl, c_incl;
    double s_decl, s_incl;

    mag->mag_scale = mag_scale;
    mag->mag_decl = mag_decl;
    mag->mag_incl = mag_incl;
    mag->mag_noise_std_dev = mag_noise_std_dev;

    c_decl = cos(deg2rad(mag->mag_decl));
    c_incl = cos(deg2rad(mag->mag_incl));
    s_decl = sin(deg2rad(mag->mag_decl));
    s_incl = sin(deg2rad(mag->mag_incl));

    mag->mag_field[0] = mag->mag_scale*c_incl*c_decl;
    mag->mag_field[1] = mag->mag_scale*c_incl*s_decl;
    mag->mag_field[2] = mag->mag_scale*s_incl;
}

void update_mag(MagSensor *mag, double dcm_be[3][3])
{
    double c_decl, c_incl;
    double s_decl, s_incl;

    double mag_e[3];

    c_decl = cos(deg2rad(mag->mag_decl));
    c_incl = cos(deg2rad(mag->mag_incl));
    s_decl = sin(deg2rad(mag->mag_decl));
    s_incl = sin(deg2rad(mag->mag_incl));

    mag_e[0] = c_incl*c_decl;
    mag_e[1] = c_incl*s_decl;
    mag_e[2] = s_incl;

    earth_to_body_rotation(dcm_be, mag_e, mag->mag_field);

    mag->mag_field[0] = mag->mag_scale*mag->mag_field[0];// + zero_mean_noise(mag->mag_noise_std_dev);
    mag->mag_field[1] = mag->mag_scale*mag->mag_field[1];// + zero_mean_noise(mag->mag_noise_std_dev);
    mag->mag_field[2] = mag->mag_scale*mag->mag_field[2];// + zero_mean_noise(mag->mag_noise_std_dev);
}

void init_baro(BaroSensor *baro, double temperature)
{
    baro->temperature = temperature;
    baro->pressure_alt = HOME_ALT;
    baro->pressure = Pb * pow((Tb / (Tb + (Lb * HOME_ALT))), ((GRAVITY * M) / (R * Lb)));
}

void update_baro(BaroSensor *baro, double noisy_alt)
{
    baro->pressure_alt = noisy_alt;
    baro->pressure = Pb * pow((Tb / (Tb + (Lb * noisy_alt))), ((GRAVITY * M) / (R * Lb)));
}