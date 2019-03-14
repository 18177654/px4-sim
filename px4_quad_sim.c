#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>

#include "quad_parameters.h"
#include "utils.h"
#include "px4_sim_communication.h"
#include "px4_quad_sim.h"

// Function prototypes.
void testQuadDynamics();
void printStates(Quad *quad);
void printSensors(Quad *quad);
void sighandler(int sig);

// Function definitions.
int init_sim()
{
    int i;
    double inertia[3];
    double c_d[3];

    // Initialise quad model
    inertia[0] = I_XX;
    inertia[1] = I_YY;
    inertia[2] = I_ZZ;

    for(i = 0 ; i < 3 ; i++)
        c_d[i] = C_D;

    init_quad(&quad, MASS, inertia, D_ARM, R_D, c_d, THRUST_TC, THRUST_HOVER_NORM * THRUST_MAX_FORCE, HOME_YAW);
    init_quad_sensors(&quad, GPS_EPH, GPS_EPV, GPS_FIX, GPS_NUM_SATS, GPS_LAT_LON_NOISE, GPS_ALT_NOISE, GPS_SPEED_NOISE, IMU_ACC_NOISE, IMU_GYRO_NOISE, MAG_NOISE, BARO_NOISE);

    // Initialise PX4 sim
    return init_px4_sim(SENSOR_FREQ, GPS_FREQ, HIL, THRUST_HOVER_NORM, THRUST_MAX_FORCE);
}

int advance_sim(uint64_t time_usec, uint64_t prev_time_usec, double wind_vel_e[3])
{
    int i;
    int result;
    double latlonalt[3];
    double thrust_commands[4];
    double dt;

    // Send HIL MAVLink messages
    ned_to_latlonalt(quad.state.pos_e, latlonalt, HOME_LAT, HOME_LON, HOME_ALT);
    result = send_hil_messages(time_usec, quad.state.quat, quad.state.euler_rates, quad.state.acc_b, quad.state.dcm_be, quad.state.vel_e, latlonalt, quad.sensors.gps.lat_lon_alt, quad.sensors.gps.gps_speed, quad.sensors.gps.ground_speed, quad.sensors.gps.cog, quad.sensors.gps.eph, quad.sensors.gps.epv, quad.sensors.gps.fix, quad.sensors.gps.visible_sats, quad.sensors.imu.acc, quad.sensors.imu.gyro, quad.sensors.mag.mag_field, quad.sensors.baro.pressure, quad.sensors.baro.diff_pressure, quad.sensors.baro.pressure_alt, quad.sensors.baro.temperature);
    if(result < 0)
        return result;

    if(result > 0)
    {
        dt = (double)((time_usec - prev_time_usec) * 1e-6);
        
        // Gaurd against large dt's
        // if(dt >= 2 * (double)(1.0 / SENSOR_FREQ))
        // {
        //     dt = 2.0 * (double)(1.0 / SENSOR_FREQ);
        //     LOG_MSG("changed dt: %2.4f ms\n", dt);
        // }

        get_thrust_commands_force(thrust_commands);
        update_quad(&quad, thrust_commands, wind_vel_e, dt);
    }

    return result;
}

int main()
{
    int i;
    uint64_t cur_time;
    uint64_t prev_time;
    double wind_vel_e[3];
    int result;
    double excess_time;

    signal(SIGINT, sighandler); // Exit simulation if Ctrl-C is pressed

    // Connect to PX4 and initialize physics sim
    LOG_MSG("Connecting to PX4...\n");
    result = init_sim();
    if(result < 0)
    {
        LOG_MSG("Failed! Result: %d\n", result);
        return 0;
    }
    LOG_MSG("Connected!\n");

    // Initialise variables
    if(HIL)
        cur_time = get_time_usec(); // in HIL, use real time.
    else
        cur_time = 0; // in SIL, use simulation time
    prev_time = cur_time;

    for(i = 0 ; i < 3 ; i++)
        wind_vel_e[i] = 0.0f;

    // Loop
    while(result >= 0)
    {
        // Update the physics sim with a fixed time
        if(HIL)
            cur_time = get_time_usec();
        else
            cur_time += 1000000.0/SENSOR_FREQ;
        
        result = advance_sim(cur_time, prev_time, wind_vel_e);
        if(result > 0)
            prev_time = cur_time;
        
        // Poll for MAVLink messages: receive actuator controls from PX4 and when in HIL - send messages from autopilot to GCS and vice-versa.
        result = pollMavlinkMessage();
        if(result < 0)
            return result;
        
        // Sleep
        if(HIL)
            usleep(10); // 10 us
        else
            usleep((1000000.0/SENSOR_FREQ)/SIL_SPEED_FACTOR);
    }

    LOG_MSG("Simulation terminated...\n");
    return 0;
}

void sighandler(int sig)
{
   LOG_MSG("\nExit simulation...\n");
   disconnect_sim();
   exit(1);
}

// Only for testing:
void testQuadDynamics()
{
    int i;
    double dt = 0;
    double sim_time = 0;
    double stop_time = 1.0;//(double)(1.0 / SENSOR_FREQ);
    double wind_vel_e[3];
    double inertia[3];
    double c_d[3];
    double thrust_commands[4];
    double forces[3];
    double moments[3];

    // Initialise quad model
    inertia[0] = I_XX;
    inertia[1] = I_YY;
    inertia[2] = I_ZZ;

    for(i = 0 ; i < 3 ; i++)
        c_d[i] = C_D;

    init_quad(&quad, MASS, inertia, D_ARM, R_D, c_d, THRUST_TC, THRUST_HOVER_NORM * THRUST_MAX_FORCE, HOME_YAW);
    init_quad_sensors(&quad, GPS_EPH, GPS_EPV, GPS_FIX, GPS_NUM_SATS, GPS_LAT_LON_NOISE, GPS_ALT_NOISE, GPS_SPEED_NOISE, IMU_ACC_NOISE, IMU_GYRO_NOISE, MAG_NOISE, BARO_NOISE);

    // Initialise wind and thrust
    // for(i = 0 ; i < 3 ; i++)
    //     wind_vel_e[i] = 0.0f;
    
    // thrust_commands[0] = 0.3 * (0.5 * MASS * GRAVITY);
    // thrust_commands[1] = 0.6 * (0.5 * MASS * GRAVITY);
    // thrust_commands[2] = 0.7 * (0.5 * MASS * GRAVITY);
    // thrust_commands[3] = 0.6 * (0.5 * MASS * GRAVITY);

    // Initialise forces and moments
    forces[0] = 1;
    forces[1] = 0;
    forces[2] = 3;

    moments[0] = 0;
    moments[1] = 2;
    moments[2] = 0;

    while(sim_time <= stop_time)
    {
        // update_quad(&quad, thrust_commands, wind_vel_e, dt);
        // printStates(&quad);
        // printSensors(&quad);

        six_dof(dt, &quad, forces, moments);
        printStates(&quad);

        dt = (double)(1.0 / SENSOR_FREQ);
        sim_time += dt;
    }
}

void printStates(Quad *quad)
{
    int i;

    LOG_MSG("Acceleration (Body):\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.9f\t", quad->state.acc_b[i]);
    LOG_MSG("\n");

    LOG_MSG("Velocity (Body):\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.9f\t", quad->state.vel_b[i]);
    LOG_MSG("\n");

    LOG_MSG("Omega (Body):\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.9f\t", quad->state.omega_b[i]);
    LOG_MSG("\n");

    LOG_MSG("Euler:\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.9f\t", quad->state.euler[i] * (180.0 / M_PI));
    LOG_MSG("\n");

    LOG_MSG("Velocity (earth):\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.9f\t", quad->state.vel_e[i]);
    LOG_MSG("\n");

    LOG_MSG("Position (earth):\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.9f\t", quad->state.pos_e[i]);
    LOG_MSG("\n");

    LOG_MSG("\n");
}

void printSensors(Quad *quad)
{
    int i;

    LOG_MSG("Accelerometer:\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.6f\t", quad->sensors.imu.acc[i]);
    LOG_MSG("\n");

    LOG_MSG("Gyro:\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.6f\t", quad->sensors.imu.gyro[i]);
    LOG_MSG("\n");

    LOG_MSG("Barometer:\n");
    LOG_MSG("%4.6f\t", quad->sensors.baro.pressure);
    LOG_MSG("\n");

    LOG_MSG("Magnetometer:\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.6f\t", quad->sensors.mag.mag_field[i]);
    LOG_MSG("\n");

    LOG_MSG("GPS:\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.6f\t", quad->sensors.gps.lat_lon_alt[i]);
    LOG_MSG("\n");

    LOG_MSG("GPS Speed:\n");
    for(i = 0 ; i < 3 ; i++)
        LOG_MSG("%4.6f\t", quad->sensors.gps.gps_speed[i]);
    LOG_MSG("\n\n");
}