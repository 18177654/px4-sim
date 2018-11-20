#ifndef PX4_SIM_COMMUNICATION_H_
#define PX4_SIM_COMMUNICATION_H_

#include <inttypes.h>
#include "mavlink/common/mavlink.h"

// Constants.
#define ERROR_CONNECT -1
#define ERROR_TTY -2
#define ERROR_GET_CONFIG -3
#define ERROR_APPLY_CONFIG -4
#define ERROR_SERVER -5
#define ERROR_RECV_TIMEOUT -6
#define ERROR_SEND_TIMEOUT -7
#define ERROR_BIND -8

#define ERROR_READ_FAIL -9
#define ERROR_READ_TIMEOUT -10

#define ERROR_WRITE_FAIL -11

// Global variables.
mavlink_message_t message; // Last read message

// Function prototypes.
int init_px4_sim(uint16_t sensor_freq, uint16_t gps_freq, int hil_enabled, float thrust_hover_norm, float max_thrust_force);
void disconnect_sim();

int send_hil_messages(uint64_t time_usec, double q[4], double lat_lon_alt[3], double vel_e[3], double vel, double cog, double eph, double epv, int fix_type, int num_sats, double acc_b[3], double gyro[3], double mag[3], double pressure, double temperature);

int pollMavlinkMessage();

void get_thrust_commands_force(double thrust_commands[4]);

#endif