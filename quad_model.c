#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "quad_model.h"

// Private globals.
Quad *_tempQuad;
double *_tempForces;
double *_tempMoments;
double _tempThrustCommand;

// Private function prototypes
void update_sensors(Quad *quad);

double u_dot(double uvwpqr_euler_ned[], int len);
double v_dot(double uvwpqr_euler_ned[], int len);
double w_dot(double uvwpqr_euler_ned[], int len);
double p_dot(double uvwpqr_euler_ned[], int len);
double q_dot(double uvwpqr_euler_ned[], int len);
double r_dot(double uvwpqr_euler_ned[], int len);

double phi_dot(double uvwpqr_euler_ned[], int len);
double theta_dot(double uvwpqr_euler_ned[], int len);
double psi_dot(double uvwpqr_euler_ned[], int len);
double n_dot(double uvwpqr_euler_ned[], int len);
double e_dot(double uvwpqr_euler_ned[], int len);
double d_dot(double uvwpqr_euler_ned[], int len);

double thrust_dot(double thrust[], int len);

// Public function definitions
void init_quad(Quad *quad, double mass, double inertia[3], double d, double r_d, double r_ld, double c_d[3], double thrust_tc, double throttle_hover)
{
    int i;
    int j;

    quad->mass = mass;
    quad->d = d;
    quad->r_d = r_d;
    quad->r_ld = r_ld;
    quad->thrust_tc = thrust_tc;
    
    for(i = 0 ; i < 4 ; i++)
        quad->thrust[i] = throttle_hover / 4;

    for(i = 0 ; i < 3 ; i++)
    {
        quad->inertia[i] = inertia[i];
        quad->c_d[i] = c_d[i];

        quad->state.acc_e[i] = 0.0f;
        quad->state.acc_b[i] = 0.0f;
        quad->state.vel_b[i] = 0.0f;
        quad->state.vel_e[i] = 0.0f;
        quad->state.pos_e[i] = 0.0f;
        quad->state.alpha_b[i] = 0.0f;
        quad->state.omega_b[i] = 0.0f;
        quad->state.euler_rates[i] = 0.0f;
        quad->state.euler[i] = 0.0f;

        for(j = 0 ; j < 3 ; j++)
        {
            if(i == j)
            {
                quad->state.dcm_be[i][j] = 1.0f;
                continue;
            }
            quad->state.dcm_be[i][j] = 0.0f;
        }
    }
}

void init_quad_sensors(Quad *quad, double eph, double epv, double fix, double visible_sats, double lat_lon_noise_std_dev, double alt_noise_std_dev, double speed_noise_std_dev, double acc_noise_std_dev, double gyro_noise_std_dev, double mag_decl, double mag_incl, double mag_scale, double mag_noise_std_dev, double temperature)
{
    init_gps(&(quad->sensors.gps), eph, epv, fix, visible_sats, lat_lon_noise_std_dev, alt_noise_std_dev, speed_noise_std_dev);
    init_imu(&(quad->sensors.imu), acc_noise_std_dev, gyro_noise_std_dev);
    init_mag(&(quad->sensors.mag), mag_decl, mag_incl, mag_scale, mag_noise_std_dev);
    init_baro(&(quad->sensors.baro), temperature);
}

void six_dof(double dt, Quad *quad, double forces[3], double moments[3])
{
    int i;
    int len = 12;
    double uvwpqr_euler_ned[len];
    double integ_uvwpqr_euler_ned[len];
    functiontype uvwpqr_euler_ned_dot[] = {u_dot, v_dot, w_dot, p_dot, q_dot, r_dot, phi_dot, theta_dot, psi_dot, n_dot, e_dot, d_dot};
    QuadState *quad_state;

    quad_state = &(quad->state);

    // Store quad, forces and moments to use in 6DOF functions.
    _tempQuad = quad;
    _tempForces = forces;
    _tempMoments = moments;

    // Store the current states.
    for(i = 0 ; i < 3 ; i++)
    {
        uvwpqr_euler_ned[i] = quad_state->vel_b[i];
        uvwpqr_euler_ned[i+3] = quad_state->omega_b[i];
        uvwpqr_euler_ned[i+6] = quad_state->euler[i];
        uvwpqr_euler_ned[i+9] = quad_state->pos_e[i];
    }

    // Integrate states using Runge Kutta.
    integrate_rk4(uvwpqr_euler_ned_dot, uvwpqr_euler_ned, integ_uvwpqr_euler_ned, len, dt);
    for(i = 0 ; i < 3 ; i++) 
    {
        quad_state->vel_b[i] = integ_uvwpqr_euler_ned[i];
        quad_state->omega_b[i] = integ_uvwpqr_euler_ned[i+3];
        quad_state->euler[i] = wrap_angle_pi(integ_uvwpqr_euler_ned[i+6]);
        quad_state->pos_e[i] = integ_uvwpqr_euler_ned[i+9];
    }

    // Calculate and store euler rates.
    quad_state->euler_rates[0] = phi_dot(integ_uvwpqr_euler_ned, len);
    quad_state->euler_rates[1] = theta_dot(integ_uvwpqr_euler_ned, len);
    quad_state->euler_rates[2] = psi_dot(integ_uvwpqr_euler_ned, len);

    // Calculate and store body angular acceleration.
    quad_state->alpha_b[0] = p_dot(integ_uvwpqr_euler_ned, len);
    quad_state->alpha_b[1] = q_dot(integ_uvwpqr_euler_ned, len);
    quad_state->alpha_b[2] = r_dot(integ_uvwpqr_euler_ned, len);

    // Calculate and store the DCM and rotate body velocity to earth velocity.
    calc_dcm_be(quad_state->euler, quad_state->dcm_be);
    body_to_earth_rotation(quad_state->dcm_be, quad_state->vel_b, quad_state->vel_e);

    // Calculate and store the acceleration.
    quad_state->acc_b[0] = u_dot(integ_uvwpqr_euler_ned, len);
    quad_state->acc_b[1] = v_dot(integ_uvwpqr_euler_ned, len);
    quad_state->acc_b[2] = w_dot(integ_uvwpqr_euler_ned, len);
    body_to_earth_rotation(quad_state->dcm_be, quad_state->acc_b, quad_state->acc_e);
}

void forces_moments_aerodynamic_model(Quad *quad, double wind_vel_e[3], double forces[], double moments[])
{
    int i;
    double scaled_wind[3];
    double wind_vel_b[3];
    double rel_vel[3];

    // Scale wind velocity in inertial frame with altitude.
    for(i = 0 ; i < 3 ; i++)
        scaled_wind[i] = pow(abs(quad->state.pos_e[2])/10.0, ETA_W) * wind_vel_e[i];

    // Convert scaled wind to body axes.
    earth_to_body_rotation(quad->state.dcm_be, scaled_wind, wind_vel_b);

    // Calculate relative velocity.
    for(i = 0 ; i < 3 ; i++)
        rel_vel[i] = -quad->state.vel_b[i] + wind_vel_b[i];

    // Calculate aerodynamic forces.
    for(i = 0 ; i < 3 ; i++)
        forces[i] = 0.5 * RHO * abs(rel_vel[i]) * rel_vel[i] * quad->c_d[i];

    for(i = 0 ; i < 3 ; i++)
        moments[i] = 0;
}

void forces_moments_thrust_model(double dt, double thrust_commands[4], Quad *quad, double forces[], double moments[])
{
    int i;
    double integ_thrust;
    functiontype thrust_dot_f[] = {thrust_dot};

    _tempQuad = quad;

    for(i = 0 ; i < 4 ; i++)
    {
        _tempThrustCommand = thrust_commands[i];
        integrate_rk4(thrust_dot_f, &(quad->thrust[i]), &integ_thrust, 1, dt);
        quad->thrust[i] = integ_thrust;
    }

    forces[0] = 0.0f;
    forces[1] = 0.0f;
    forces[2] = -(quad->thrust[0]+quad->thrust[1]+quad->thrust[2]+quad->thrust[3]);

    moments[0] = quad->d * (quad->thrust[3] - quad->thrust[1]);
    moments[1] = quad->d * (quad->thrust[0] - quad->thrust[2]);
    moments[2] = quad->r_d * (-quad->thrust[0]+quad->thrust[1]-quad->thrust[2]+quad->thrust[3]) / quad->r_ld;
}

void forces_moments_gravity_model(Quad *quad, double forces[], double moments[])
{
    int i;

    for(i = 0 ; i < 3 ; i++)
    {
        forces[i] = quad->mass * GRAVITY * quad->state.dcm_be[i][2];
        moments[i] = 0;
    }
}

void get_att_quaternion(Quad *quad, double q[4])
{
    q[0] = 0.0f;
    q[1] = 0.0f;
    q[2] = 0.0f;
    q[3] = 0.0f;
}

void update_quad(Quad *quad, double thrust_commands[4], double wind_vel_e[3], double dt)
{
    int i;
    double forces[3], forces_aero[3], forces_thrust[3], forces_gravity[3];
    double moments[3], moments_aero[3], moments_thrust[3], moments_gravity[3];

    forces_moments_aerodynamic_model(quad, wind_vel_e, forces_aero, moments_aero);
    forces_moments_gravity_model(quad, forces_gravity, moments_gravity);
    forces_moments_thrust_model(dt, thrust_commands, quad, forces_thrust, moments_thrust);

    for(i = 0 ; i < 3 ; i++)
    {
        forces[i] = forces_thrust[i] + forces_gravity[i] + forces_aero[i];
        moments[i] = moments_thrust[i] + moments_gravity[i] + moments_aero[i];
    }

    six_dof(dt, quad, forces, moments);
    update_sensors(quad);
}

void update_sensors(Quad *quad)
{
    update_gps(&(quad->sensors.gps), quad->state.pos_e, quad->state.vel_e);
    update_imu(&(quad->sensors.imu), quad->state.acc_b, quad->state.omega_b, quad->state.dcm_be);
    update_mag(&(quad->sensors.mag), quad->state.dcm_be);
    update_baro(&(quad->sensors.baro), quad->sensors.gps.lat_lon_alt[2]);
}

// Differential functions.
double u_dot(double uvwpqr_euler_ned[], int len)
{
    return _tempForces[0]/_tempQuad->mass + uvwpqr_euler_ned[1]*uvwpqr_euler_ned[5] - uvwpqr_euler_ned[2]*uvwpqr_euler_ned[4];
}

double v_dot(double uvwpqr_euler_ned[], int len)
{
    return _tempForces[1]/_tempQuad->mass - uvwpqr_euler_ned[0]*uvwpqr_euler_ned[5] + uvwpqr_euler_ned[2]*uvwpqr_euler_ned[3];
}

double w_dot(double uvwpqr_euler_ned[], int len)
{
    return _tempForces[2]/_tempQuad->mass + uvwpqr_euler_ned[0]*uvwpqr_euler_ned[4] - uvwpqr_euler_ned[1]*uvwpqr_euler_ned[3];
}

double p_dot(double uvwpqr_euler_ned[], int len)
{
    return (_tempMoments[0] - uvwpqr_euler_ned[4]*uvwpqr_euler_ned[5]*(_tempQuad->inertia[2] - _tempQuad->inertia[1]))/_tempQuad->inertia[0];
}

double q_dot(double uvwpqr_euler_ned[], int len)
{
    return (_tempMoments[1] - uvwpqr_euler_ned[3]*uvwpqr_euler_ned[5]*(_tempQuad->inertia[0] - _tempQuad->inertia[2]))/_tempQuad->inertia[1];
}

double r_dot(double uvwpqr_euler_ned[], int len)
{
    return (_tempMoments[2] - uvwpqr_euler_ned[3]*uvwpqr_euler_ned[4]*(_tempQuad->inertia[1] - _tempQuad->inertia[0]))/_tempQuad->inertia[2];
}

double phi_dot(double uvwpqr_euler_ned[], int len)
{
    return uvwpqr_euler_ned[3] + sin(uvwpqr_euler_ned[6])*tan(uvwpqr_euler_ned[7])*uvwpqr_euler_ned[4] + cos(uvwpqr_euler_ned[6])*tan(uvwpqr_euler_ned[7])*uvwpqr_euler_ned[5];
}

double theta_dot(double uvwpqr_euler_ned[], int len)
{
    return cos(uvwpqr_euler_ned[6])*uvwpqr_euler_ned[4] - sin(uvwpqr_euler_ned[7])*uvwpqr_euler_ned[5];
}

double psi_dot(double uvwpqr_euler_ned[], int len)
{
    return (sin(uvwpqr_euler_ned[6])/cos(uvwpqr_euler_ned[7]))*uvwpqr_euler_ned[4] + (cos(uvwpqr_euler_ned[6])/cos(uvwpqr_euler_ned[7]))*uvwpqr_euler_ned[5];
}

double n_dot(double uvwpqr_euler_ned[], int len)
{
    double c_phi, c_theta, c_psi;
    double s_phi, s_theta, s_psi;

    c_phi = cos(uvwpqr_euler_ned[6]);
    c_theta = cos(uvwpqr_euler_ned[7]);
    c_psi = cos(uvwpqr_euler_ned[8]);

    s_phi = sin(uvwpqr_euler_ned[6]);
    s_theta = sin(uvwpqr_euler_ned[7]);
    s_psi = sin(uvwpqr_euler_ned[8]);

    return c_psi*c_theta*uvwpqr_euler_ned[0] + (c_psi*s_theta*s_phi - s_psi*c_phi)*uvwpqr_euler_ned[1] + (c_psi*s_theta*c_phi + s_psi*s_phi)*uvwpqr_euler_ned[2];
}

double e_dot(double uvwpqr_euler_ned[], int len)
{
    double c_phi, c_theta, c_psi;
    double s_phi, s_theta, s_psi;

    c_phi = cos(uvwpqr_euler_ned[6]);
    c_theta = cos(uvwpqr_euler_ned[7]);
    c_psi = cos(uvwpqr_euler_ned[8]);

    s_phi = sin(uvwpqr_euler_ned[6]);
    s_theta = sin(uvwpqr_euler_ned[7]);
    s_psi = sin(uvwpqr_euler_ned[8]);

    return s_psi*c_theta*uvwpqr_euler_ned[0] + (s_psi*s_theta*s_phi + c_psi*c_phi)*uvwpqr_euler_ned[1] + (s_psi*s_theta*c_phi - c_psi*s_phi)*uvwpqr_euler_ned[2];
}

double d_dot(double uvwpqr_euler_ned[], int len)
{
    double c_phi, c_theta, c_psi;
    double s_phi, s_theta, s_psi;

    c_phi = cos(uvwpqr_euler_ned[6]);
    c_theta = cos(uvwpqr_euler_ned[7]);
    c_psi = cos(uvwpqr_euler_ned[8]);

    s_phi = sin(uvwpqr_euler_ned[6]);
    s_theta = sin(uvwpqr_euler_ned[7]);
    s_psi = sin(uvwpqr_euler_ned[8]);

    return -s_theta*uvwpqr_euler_ned[0] + s_psi*c_theta*uvwpqr_euler_ned[1] + (c_theta*c_phi)*uvwpqr_euler_ned[2];
}

double thrust_dot(double thrust[], int len)
{
    return (-thrust[0] + _tempThrustCommand)/_tempQuad->thrust_tc;
}