#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "utils.h"

double deg2rad(double deg)
{
    return deg * M_PI / 180.0f;
}

double rad2deg(double rad)
{
    return rad * 180.0f / M_PI;
}

double zero_mean_noise(double std_dev)
{
    return rand_gauss(0.0, std_dev);
}

// Gauss distribution - polar method
double rand_gauss(double mean, double std_dev)
{
    static int first = 1;
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if(first)
    {
        srand(time(NULL));
        first = 0;
    }

    if (call == 1)
    {
        call = !call;
        return (mean + std_dev * (double)X2);
    }

    do
    {
        U1 = -1 + ((double)rand() / RAND_MAX) * 2;
        U2 = -1 + ((double)rand() / RAND_MAX) * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mean + std_dev * (double)X1);
}

void calc_dcm_be(double euler[3], double dcm_be[3][3])
{
    double c_phi, c_theta, c_psi;
    double s_phi, s_theta, s_psi;

    c_phi = cos(euler[0]);
    c_theta = cos(euler[1]);
    c_psi = cos(euler[2]);

    s_phi = sin(euler[0]);
    s_theta = sin(euler[1]);
    s_psi = sin(euler[2]);

    dcm_be[0][0] = c_psi*c_theta;
    dcm_be[0][1] = s_psi*c_theta;
    dcm_be[0][2] = -s_theta;
    dcm_be[1][0] = c_psi*s_theta*s_phi - s_psi*c_phi;
    dcm_be[1][1] = s_psi*s_theta*s_phi + c_psi*c_phi;
    dcm_be[1][2] = c_theta*s_phi;
    dcm_be[2][0] = c_psi*s_theta*c_phi + s_psi*s_phi;
    dcm_be[2][1] = s_psi*s_theta*c_phi - c_psi*s_phi;
    dcm_be[2][2] = c_theta*c_phi;
}

void body_to_earth_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3])
{
    int i, j;

    for(i = 0 ; i < 3 ; i++)
    {   
        rotate_to[i] = 0;
        for(j = 0 ; j < 3 ; j++)
        {
            rotate_to[i] += (dcm_be[j][i] * rotate_from[j]);
        }
    }
}

void earth_to_body_rotation(double dcm_be[3][3], double rotate_from[3], double rotate_to[3])
{
    int i, j;

    for(i = 0 ; i < 3 ; i++)
    {   
        rotate_to[i] = 0;
        for(j = 0 ; j < 3 ; j++)
        {
            rotate_to[i] += (dcm_be[i][j] * rotate_from[j]);
        }
    }
}

double wrap_angle_2pi(double angle)
{
    double wrapped;

    if(angle >= 0 && angle <= 2*M_PI)
        return angle;

    // Wrap between 0 and 2pi.
    wrapped = angle - (floor(angle/(2*M_PI)))*2*M_PI;

    return wrapped;
}

double wrap_angle_pi(double angle)
{
    double wrapped;

    if(angle >= -M_PI && angle <= M_PI)
        return angle;

    // Wrap between 0 and 2pi.
    wrapped = wrap_angle_2pi(angle);

    // Make sure wrapped angle is between pi and -pi
    if(wrapped > M_PI)
        wrapped = wrapped - (2.0*M_PI);

    return wrapped;
}

void integrate_rk4(functiontype derivs[], double states[], double integ[], int len, double h)
{
    int i;
    double k1[len], k2[len], k3[len], k4[len];
    double states_k1[len], states_k2[len], states_k3[len];

    // k1
    for(i = 0 ; i < len ; i++)
    {
        k1[i] = h*derivs[i](states, len);
        states_k1[i] = states[i] + k1[i]/2.0;
    }

    // k2
    for(i = 0 ; i < len ; i++)
    {
        k2[i] = h*derivs[i](states_k1, len);
        states_k2[i] = states[i] + k2[i]/2.0;
    }

    // k3
    for(i = 0 ; i < len ; i++)
    {
        k3[i] = h*derivs[i](states_k2, len);
        states_k3[i] = states[i] + k3[i];
    }

    // k4
    for(i = 0 ; i < len ; i++)
        k4[i] = h*derivs[i](states_k3, len);

    // weigh
    for(i = 0 ; i < len ; i++)
        integ[i] = states[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
}

double integrate_euler(double sum, double val, double dt)
{
    return sum + val*dt;
}

uint64_t get_time_usec()
{
    static struct timeval _time_stamp;
    gettimeofday(&_time_stamp, NULL);
    return _time_stamp.tv_sec * 1000000 + _time_stamp.tv_usec;
}