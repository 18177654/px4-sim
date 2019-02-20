#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <termios.h>
#include <arpa/inet.h>
#include <netinet/tcp.h>
#include <sys/socket.h>
#include <poll.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#include "px4_sim_communication.h"
#include "utils.h"

// Constants.
#define ENABLED 1
#define DISABLED 0

#define DEBUG_SIM ENABLED
// #define DEBUG_SIM DISABLED

#define DEFAULT_HIL_COM_PORT "/dev/ttyACM0"
#define DEFAULT_HIL_BAUD B921600

#define DEFAULT_HIL_GCS_SERVER "127.0.0.1"
#define DEFAULT_HIL_GCS_PORT 14550

#define DEFAULT_SIL_SERVER "127.0.0.1"
#define DEFAULT_SIL_PORT 4560

#define HIL_MSG_TIME_ERROR 0 // us

#define TIMEOUT 1000 // ms
#define MSG_ID_TIMEOUT 2000 // in ms, the timeout when waiting for a specific message

#define POLL_TIMEOUT 0 // ms

#define HIL_POLL_AP 0 // autopilot poll index
#define HIL_POLL_GCS 1 // ground control station poll index

#define SIL_POLL_AP 2 // autopilot poll index

// Global variables.
int hil; // 0 = SIL, 1 = HIL

char hil_com_port[15];
int hil_baud;

char hil_gcs_server[15];
int hil_gcs_port;

char sil_server[15];
int sil_port;

int fc_serial_fd; // serial connection to flight controller
struct termios fc_serial_config; // serial configuration

int gcs_udp_fd; // udp connection to GCS
struct sockaddr_in gcs_udp_addr_in; // udp address configuration

int sitl_tcp_fd; // tcp connection to PX4 SITL
struct sockaddr_in sitl_tcp_addr_in; // tcp address configuration

mavlink_hil_actuator_controls_t hil_actuators; // Last HIL actuators message received.

struct pollfd poll_fds[3];

float hil_state_freq; // Frequency at which HIL_STATE should be sent

uint64_t time_usec_start; // The time when the simulation started
uint64_t time_usec; // The current simulation time in us

uint32_t hil_sensor_startup_delay;
uint32_t hil_gps_startup_delay;

uint16_t hil_sensor_freq;
uint16_t hil_gps_freq;

int8_t mav_sys_id; // The system ID of the autopilot
int8_t mav_comp_id; // The component ID of the simulator

bool mav_inited; // Indicates whether the MAVLink communication is initialised, i.e when the sim has received a heartbeat from the autopilot to obtain its system ID
bool armed; // Indicates whether the autopilot is armed
bool hover_reached; // Indicates whether the hover throttle has been reached or not.

float thrust_hover;
float max_thrust;

// Function prototypes.
void init_globals(int hil_enabled, float thrust_hover_norm, float max_thrust_force);

int connect_sim();

int connect_serial(int *fc_serial_fd, struct termios *fc_serial_config, char *com_port, int baud);
int configure_serial(int *fc_serial_fd, struct termios *fc_serial_config, int baud);

int connect_ip(int *ip_fd, struct sockaddr_in *ip_addr_in, char *server, int port, bool is_tcp);

int get_poll_index();

int read_bytes(uint8_t *buffer, uint16_t len, uint8_t poll_id);
int write_bytes(uint8_t *buffer, uint16_t len, uint8_t poll_id);

int read_message(mavlink_message_t *message, uint8_t poll_id);
int handleMessage(mavlink_message_t *message, int poll_index);

int init_autopilot();

void show_message(mavlink_message_t *message);
void show_heartbeat(mavlink_message_t *message);
void show_statustext(mavlink_message_t *message);
void show_command_long(mavlink_message_t *message);
void show_actuators(mavlink_message_t *message);

int write_message(mavlink_message_t *message, int poll_id);
int send_heartbeat();
int send_cmd_ack(uint16_t cmd, uint8_t result);
int send_cmd_set_mode();
int send_cmd_arm(uint8_t arm);
int send_hil_state(uint64_t time_usec, double q[4], double euler_rates[3], double lat_lon_alt[3], double vel_e[3], double acc_b[3]);
int send_hil_sensors(uint64_t time_usec, double acc_b[3], double gyro[3], double mag[3], double pressure, double alt, double temperature);
int send_hil_gps(uint64_t time_usec, double lat_lon_alt[3], double vel_e[3], double vel, double cog, double eph, double epv, int fix_type, int num_sats);

// Function definitions.

int init_px4_sim(uint16_t sensor_freq, uint16_t gps_freq, int hil_enabled, float thrust_hover_norm, float max_thrust_force)
{
    // Initialize globals.
    init_globals(hil_enabled, thrust_hover_norm, max_thrust_force);

    // Connect to the device.
    int result = connect_sim();
    if (result < 0)
        return result;

    // Wait for autopilot heartbeat.
    result = init_autopilot();
    if(result < 0)
        return result;
}

/* Function: init_globals ===============================================
 * Abstract:
 *      Initialize all the global variables.
 * 
 * Parameters:
 *      hil_enabled - 1 for HIL; 0 for SIL
 *      thrust_hover - the thrust at which the quad hovers (0.0 - 1.0)
 */
void init_globals(int hil_enabled, float thrust_hover_norm, float max_thrust_force)
{
    time_usec_start = get_time_usec();

    hil = hil_enabled;

    if(hil)
    {
        strcpy(hil_com_port, DEFAULT_HIL_COM_PORT);
        hil_baud = DEFAULT_HIL_BAUD;

        strcpy(hil_gcs_server, DEFAULT_HIL_GCS_SERVER);
        hil_gcs_port = DEFAULT_HIL_GCS_PORT;

        hil_sensor_startup_delay = 5000;
        hil_gps_startup_delay = 5000;
    }
    else
    {
        strcpy(sil_server, DEFAULT_SIL_SERVER);
        sil_port = DEFAULT_SIL_PORT;

        hil_sensor_startup_delay = 0;
        hil_gps_startup_delay = 0;
    }

    mav_sys_id = -1;
    mav_comp_id = -1;

    hil_state_freq = -1;

    time_usec = 0;

    mav_inited = false;
    armed = false;
    hover_reached = false;

    thrust_hover = thrust_hover_norm;
    max_thrust = max_thrust_force;

    // Set all HIL actuators to hover.
    hil_actuators.time_usec = time_usec_start;
    for(uint8_t i = 0 ; i < 16 ; i++)
        hil_actuators.controls[i] = thrust_hover;
}

/* Function: connect_sim ===============================================
 * Abstract:
 *      Connects to and configures the flight controller.
 *
 * Returns:
 *      - 0: Connection and configuration successfull.
 *      - otherwise: an error occurred.
 */
int connect_sim()
{
    int result = 0;

    if(hil)
    {
        result = connect_serial(&fc_serial_fd, &fc_serial_config, hil_com_port, hil_baud);

        poll_fds[HIL_POLL_AP].fd = fc_serial_fd;
        poll_fds[HIL_POLL_AP].events = 0; // clear
        poll_fds[HIL_POLL_AP].events |= POLLIN; // poll read
        if(result >= 0)
        {
            result = connect_ip(&gcs_udp_fd, &gcs_udp_addr_in, hil_gcs_server, hil_gcs_port, false);

            poll_fds[HIL_POLL_GCS].fd = gcs_udp_fd;
            poll_fds[HIL_POLL_GCS].events = 0; // clear
            poll_fds[HIL_POLL_GCS].events |= POLLIN; // poll read
        }
    }
    else
    {
        result = connect_ip(&sitl_tcp_fd, &sitl_tcp_addr_in, sil_server, sil_port, true);

        poll_fds[SIL_POLL_AP].fd = sitl_tcp_fd;
        poll_fds[SIL_POLL_AP].events = 0; // clear
        poll_fds[SIL_POLL_AP].events |= POLLIN; // poll read
    }

    return result;
}

/* Function: disconnect_sim ===============================================
 * Abstract:
 *    Terminates the connection.
 */
void disconnect_sim()
{
    if(hil)
    {
        if (fc_serial_fd >= 0)
        {
            close(fc_serial_fd);
            if(gcs_udp_fd >= 0)
                close(gcs_udp_fd);
        }
    }
    else
    {
        if(sitl_tcp_fd >= 0)
            close(sitl_tcp_fd);
    }
}

/* Function: connect_serial ===============================================
 * Abstract:
 *      Connects to and configures the given serial connection.
 * 
 * Parameters:
 *      - fc_serial_fd: the serial file descriptor
 *      - fc_serial_config: the serial configuarion
 *      - com_port: the COM port to connect to
 *      - baud: the baud rate at which to communicate
 * 
 * Returns:
 *      0 if no error occurred.
 */
int connect_serial(int *fc_serial_fd, struct termios *fc_serial_config, char *com_port, int baud)
{
    *fc_serial_fd = open(com_port, O_RDWR | O_NOCTTY | O_NDELAY);
    if (*fc_serial_fd == -1)
    {
        return ERROR_CONNECT;
    }

    // Finalize
    fcntl(*fc_serial_fd, F_SETFL, 0);

    // Configure port.
    return configure_serial(fc_serial_fd, fc_serial_config, baud);
}

/* Function: configure_serial ===============================================
 * Abstract:
 *      Configures the given serial connection.
 * 
 * Parameters:
 *      - fc_serial_fd: the serial file descriptor
 *      - fc_serial_config: the serial configuarion
 *      - baud: the baud rate at which to communicate
 * 
 * Returns:
 *      0 if no error occurred.
 */
int configure_serial(int *fc_serial_fd, struct termios *fc_serial_config, int baud)
{
    if (!isatty(*fc_serial_fd))
    {
        return ERROR_TTY;
    }
    if (tcgetattr(*fc_serial_fd, fc_serial_config) < 0)
    {
        return ERROR_GET_CONFIG;
    }

    // Input flags - Turn off input processing
    // convert break to null byte, no CR to NL translation,
    // no NL to CR translation, don't mark parity errors or breaks
    // no input parity check, don't strip high bit off,
    // no XON/XOFF software flow control
    fc_serial_config->c_iflag &= ~(IGNBRK | BRKINT | ICRNL | INLCR | PARMRK | INPCK | ISTRIP | IXON | IXOFF | IXANY);

    // Output flags - Turn off output processing
    // no CR to NL translation, no NL to CR-NL translation,
    // no NL to CR translation, no column 0 CR suppression,
    // no Ctrl-D suppression, no fill characters, no case mapping,
    // no local output processing
    fc_serial_config->c_oflag &= ~(OCRNL | ONLCR | ONLRET | ONOCR | OFILL | OPOST);

    // No line processing:
    // echo off, echo newline off, canonical mode off,
    // extended input processing off, signal chars off
    fc_serial_config->c_lflag &= ~(ECHO | ECHONL | ICANON | IEXTEN | ISIG);

    // Turn off character processing
    // clear current char size mask, no parity checking, 1 stop bit, no hardware based flow control
    // no output processing, force 8 bit input, turn on read
    fc_serial_config->c_cflag &= ~(CSIZE | PARENB | CSTOPB | CRTSCTS);
    fc_serial_config->c_cflag |= (CS8 | CREAD | CLOCAL);

    // One input byte is enough to return from read()
    // Inter-character timer off
    fc_serial_config->c_cc[VMIN] = 1;
    fc_serial_config->c_cc[VTIME] = 0;

    // Apply baudrate
    cfsetospeed(fc_serial_config, baud);
    cfsetispeed(fc_serial_config, baud);

    // Apply configuration.
    if (tcsetattr(*fc_serial_fd, TCSANOW, fc_serial_config) < 0)
    {
        return ERROR_APPLY_CONFIG;
    }

    return 0;
}

/* Function: connect_ip ===============================================
 * Abstract:
 *      Connects to and configures the given IP connection.
 * 
 * Parameters:
 *      - ip_fd: the IP file descriptor
 *      - ip_addr_in: the IP socket address configuration
 *      - server: the server to connect to
 *      - port: the port to connect to
 * 
 * Returns:
 *      0 if no error occurred.
 */
int connect_ip(int *ip_fd, struct sockaddr_in *ip_addr_in, char *server, int port, bool is_tcp)
{
    int result = 0;

    if(is_tcp) {
        int _fd;
        struct sockaddr_in _myaddr;
        socklen_t _myaddr_len;

        _fd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);

        memset(&_myaddr, 0, sizeof(_myaddr));
        _myaddr.sin_family = AF_INET;
        _myaddr.sin_addr.s_addr = inet_addr(server);
        _myaddr.sin_port = htons(port);
        _myaddr_len = sizeof(_myaddr);

        int yes = 1;
        result = setsockopt(_fd, IPPROTO_TCP, TCP_NODELAY, (char *)&yes, sizeof(yes));
        if(result < 0)
            return result;

        result = bind(_fd, (struct sockaddr*) &_myaddr, _myaddr_len);
        if(result < 0)
            return result;

        result = listen(_fd, 0);
        if(result < 0)
            return result;

        ip_addr_in->sin_family = AF_INET;
        ip_addr_in->sin_addr.s_addr = inet_addr(server);
        ip_addr_in->sin_port = htons(port);
        socklen_t socklen = sizeof(*ip_addr_in);
        *ip_fd = accept(_fd, (struct sockaddr*) ip_addr_in, &socklen);
        if(*ip_fd < 0)
            return ERROR_CONNECT;

    } else {
        *ip_fd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);

        ip_addr_in->sin_family = AF_INET;
        ip_addr_in->sin_addr.s_addr = inet_addr(server);
        ip_addr_in->sin_port = htons(port);

        if(*ip_fd < 0)
            return ERROR_CONNECT;
    }

    return result;
}

/* Function: get_poll_index ===============================================
 * Abstract:
 *      Gets the AP poll index.
 * 
 * Returns:
 *      The AP poll index.
 */
int get_poll_index()
{
    int poll_index;
    
    if(hil)
        poll_index = HIL_POLL_AP;
    else
        poll_index = SIL_POLL_AP;
    
    return poll_index;
}

/* Function: read_bytes ===============================================
 * Abstract:
 *      Reads len bytes from the device and assigns to buffer.
 * 
 * Parameters:
 *      buffer: the buffer to store the read bytes in.
 *      len: the number of bytes to read.
 *      poll_id: the ID of the file descriptor.
 * 
 * Returns:
 *      This function returns the number of bytes read, which should be len.
 *      If the returned result is not len, then an error most likely occured.
 */
int read_bytes(uint8_t *buffer, uint16_t len, uint8_t poll_id)
{
    if(poll_id == HIL_POLL_AP)
    {
        return read(fc_serial_fd, buffer, len);
    }
    else if(poll_id == SIL_POLL_AP)
    {
        socklen_t addr_len = sizeof(sitl_tcp_addr_in);
        return recvfrom(sitl_tcp_fd, buffer, len, 0, (struct sockaddr *)&sitl_tcp_addr_in, &addr_len);
    }
    else
    {
        socklen_t addr_len = sizeof(gcs_udp_addr_in);
        return recvfrom(gcs_udp_fd, buffer, len, 0, (struct sockaddr *)&gcs_udp_addr_in, &addr_len);
    }
}

/* Function: write_bytes ===============================================
 * Abstract:
 *    Sends the given buffer to the device.
 * 
 * Parameters:
 *      buffer: the buffer to send.
 *      len: the number of bytes in the buffer to send.
 *      poll_id: the ID of the file descriptor.
 * 
 * Returns:
 *    This function returns the number of bytes written, which should be len.
 *    If the returned result is not len, then an error most likely occured.
 */
int write_bytes(uint8_t *buffer, uint16_t len, uint8_t poll_id)
{
    if(poll_id == HIL_POLL_AP)
    {
        int result = write(fc_serial_fd, buffer, len);

        if (result != len)
        {
            return ERROR_WRITE_FAIL;
        }

        tcdrain(fc_serial_fd);
        return 1;
    }
    else if(poll_id == SIL_POLL_AP)
    {
        int result = sendto(sitl_tcp_fd, buffer, len, 0, (struct sockaddr *)&sitl_tcp_addr_in, sizeof(sitl_tcp_addr_in));
        
        if(result < 0)
        {
            return ERROR_WRITE_FAIL;
        }

        return result;
    }
    else
    {
        int result = sendto(gcs_udp_fd, buffer, len, 0, (struct sockaddr *)&gcs_udp_addr_in, sizeof(gcs_udp_addr_in));

        if(result < 0)
        {
            return ERROR_WRITE_FAIL;
        }

        return result;
    }
}

/* Function: read_message ===============================================
 * Abstract:
 *      This function reads from the device and parses the MAVLink message one byte at a time.
 *      This function will return 1 when a full message has been parsed, 0 otherwise.
 * 
 * Parameters:
 *      message: the MAVLink message object to store the received message in.
 *      poll_id: the ID of the file descriptor.
 * 
 * Returns:
 *      If an error ocurred, -1 will be returned.
 */
int read_message(mavlink_message_t *message, uint8_t poll_id)
{
    uint16_t len = 2048;
    
    uint8_t buffer[len];
    mavlink_status_t status;
    uint8_t msg_received = 0;

    uint8_t read_result = read_bytes(buffer, len, poll_id);
    if (read_result == 0)
    {
        return 0;
    }

    if (read_result > 0)
    {
        for(uint8_t i = 0 ; i < read_result ; i++)
        {
            msg_received = mavlink_parse_char(MAVLINK_COMM_0, buffer[i], message, &status);
            if(msg_received > 0)
            {
                break;
            }
        }
        if (msg_received > 0)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    { 
        return ERROR_READ_FAIL; // error.
    }
}

/* Function: pollMavlinkMessage ===============================================
 * Abstract:
 *      This polls to see if a MAVLink message can be received.
 *      The MAVLink message will be stored in the global variable message.
 *      This function returns the id of the message.
 * 
 * Returns:
 *      This function will return a negative number if an error occurred and 0 if no message was available.
 */
int pollMavlinkMessage()
{
    int ret = 0;
    int poll_index = get_poll_index();

    int pret = poll(&poll_fds[poll_index], 1, POLL_TIMEOUT);
    if(pret > 0)
    {
        int success = read_message(&message, poll_index);
        if(success < 0)
        {
            return success; // error
        }

        if(success == 1)
        {
            ret = handleMessage(&message, poll_index);
        }
    }

    // Check for messages from GCS
    if(hil)
    {
        pret = poll(&poll_fds[HIL_POLL_GCS], 1, POLL_TIMEOUT);
        if(pret > 0)
        {
            int success = read_message(&message, HIL_POLL_GCS);
            if(success < 0)
            {
                return success; // error
            }

            if(success == 1)
            {
                ret = handleMessage(&message, HIL_POLL_GCS);
            }
        }
    }

    return ret;
}

/* Function: handleMessage ===============================================
 * Abstract:
 *      This function handles the given message and performs the appropriate actions.
 *      This function returns the id of the message.
 * 
 * Parameters:
 *      message: the MAVLink message to handle.
 *      poll_id: the ID of the file descriptor.
 * 
 * Returns:
 *      This function will return a negative number if an error occurred.
 */
int handleMessage(mavlink_message_t *message, int poll_index)
{
#if DEBUG_SIM == ENABLED
    show_message(message);
#endif

    switch(message->msgid)
    {
        case MAVLINK_MSG_ID_HEARTBEAT:
            if(!mav_inited)
            {
                mav_inited = true;
                mav_sys_id = message->sysid;
                mav_comp_id = message->compid + 1;
#if DEBUG_SIM == ENABLED
                printf("Initialised!\n");
#endif
            }
            break;
        case MAVLINK_MSG_ID_COMMAND_LONG:
            if(hil_state_freq <= 0)
            {
                mavlink_command_long_t cmd;
                mavlink_msg_command_long_decode(message, &cmd);
                if(cmd.command == MAV_CMD_SET_MESSAGE_INTERVAL && cmd.param1 == MAVLINK_MSG_ID_HIL_STATE_QUATERNION)
                {
                    hil_state_freq = cmd.param2;
#if DEBUG_SIM == ENABLED
                    printf("HIL_STATE Interval: %f\n", hil_state_freq);
#endif
                }
            }
            break;
        case MAVLINK_MSG_ID_HIL_ACTUATOR_CONTROLS:
            mavlink_msg_hil_actuator_controls_decode(message, &hil_actuators);
            if (hil_actuators.mode != 0) {
                if ((hil_actuators.mode & 128) > 0 /* armed */) {
                    armed = true;
                } else {
                    armed = false;
                }
            }
            if(armed && !hover_reached)
            {
                float throttle = 0;
                for(uint8_t i = 0 ; i < 4 ; i++)
                {
                    if(hil_actuators.controls[i] >= 0 && hil_actuators.controls[i] <= 1)
                        throttle += hil_actuators.controls[i];
                }
                throttle /= 4;
                if(throttle >= thrust_hover)
                {
                    hover_reached = true;
                }
            }
            break;
    }

    if(message->msgid != MAVLINK_MSG_ID_HIL_ACTUATOR_CONTROLS && message->msgid != MAVLINK_MSG_ID_HIL_CONTROLS && message->msgid != MAVLINK_MSG_ID_HIL_GPS && message->msgid != MAVLINK_MSG_ID_HIL_SENSOR && message->msgid != MAVLINK_MSG_ID_HIL_STATE && message->msgid != MAVLINK_MSG_ID_HIL_STATE_QUATERNION && message->msgid != MAVLINK_MSG_ID_HIL_RC_INPUTS_RAW && message->msgid != MAVLINK_MSG_ID_HIL_OPTICAL_FLOW)
    {
        // Sim is gateway between autopilot and ground control station
        if(hil)
        {
            if(poll_index == HIL_POLL_AP)
            {
                write_message(message, HIL_POLL_GCS);
            }
            else
            {
                write_message(message, HIL_POLL_AP);
            }
        }
    }
    
    return message->msgid;
}

/* Function: init_autopilot ===============================================
 * Abstract:
 *      Sends heartbeats to the autopilot until a heartbeat is received.
 * 
 * Returns:
 *      This function returns 1 if successful.
 *      This function will return a negative number if an error occurred.
 */
int init_autopilot()
{
    int success;
    int ret = 0;

    while(!mav_inited)
    {
        ret = send_heartbeat();
        if(ret < 0)
        {
            return ret; // error.
        }

        success = pollMavlinkMessage();
        if(success < 0)
        {
            return success; // error.
        }
        
        usleep(10 * 1000); // 10 ms
    }

    ret = send_cmd_set_mode();
    if(ret < 0)
    {
        return ret; // error.
    }

    return 1;
}

/* Function: show_message ===============================================
 * Abstract:
 *      Displays the content of the message.
 * 
 * Parameters:
 *      message: the MAVLink message to show.
 */
void show_message(mavlink_message_t *message)
{
    switch (message->msgid)
    {
    case MAVLINK_MSG_ID_HEARTBEAT:
        // show_heartbeat(message);
        break;
    case MAVLINK_MSG_ID_STATUSTEXT:
        show_statustext(message);
        break;
    case MAVLINK_MSG_ID_COMMAND_LONG:
        show_command_long(message);
        break;
    case MAVLINK_MSG_ID_HIL_ACTUATOR_CONTROLS:
        show_actuators(message);
    default:
        // printf("Got message with ID: %ld\n", message->msgid);
        break;
    }
}

/* Function: show_heartbeat ===============================================
 * Abstract:
 *      Displays a message indicating that a heartbeat has been received.
 * 
 * Parameters:
 *      message: the MAVLink message to show.
 */
void show_heartbeat(mavlink_message_t *message)
{
    printf("Received heartbeat (%.4f)...\n", ((float)time_usec) / 1000000.0);
}

/* Function: show_statustext ===============================================
 * Abstract:
 *      Displays the text that has been received.
 * 
 * Parameters:
 *      message: the MAVLink message to show.
 */
void show_statustext(mavlink_message_t *message)
{
    mavlink_statustext_t statustext;
    mavlink_msg_statustext_decode(message, &statustext);
    printf("FC Message: %s\n", statustext.text);
}

/* Function: show_command_long ===============================================
 * Abstract:
 *      Displays the command that has been received.
 * 
 * Parameters:
 *      message: the MAVLink message to show.
 */
void show_command_long(mavlink_message_t *message)
{
    mavlink_command_long_t cmd;
    mavlink_msg_command_long_decode(message, &cmd);
    printf("Command received(%.4f): %d\n", ((float)time_usec) / 1000000.0, cmd.command);
}

/* Function: show_actuators ===============================================
 * Abstract:
 *      Displays the HIL actuators that has been received.
 * 
 * Parameters:
 *      message: the MAVLink message to show.
 */
void show_actuators(mavlink_message_t *message)
{
    // int i;
    // mavlink_hil_actuator_controls_t hil_actuators;
    // mavlink_msg_hil_actuator_controls_decode(message, &hil_actuators);

    // if(hil_actuators.controls[0] != hil_actuators.controls[1] || hil_actuators.controls[2] != hil_actuators.controls[3] || hil_actuators.controls[0] != hil_actuators.controls[3])
    // {
    //     for(i = 0 ; i < 4 ; i++)
    //         printf("%2.3f\t", hil_actuators.controls[i]);
    //     printf("\n");
    // }
}

/* Function: write_message ===============================================
 * Abstract:
 *      This function sends the given MAVLink message.
 * 
 * Parameters:
 *      message: the MAVLink message to send.
 *      poll_id: the ID of the file descriptor.
 * 
 * Returns:
 *      This function will return 1 when the message has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int write_message(mavlink_message_t *message, int poll_id)
{
    uint8_t buffer[2048];
    uint16_t len = mavlink_msg_to_send_buffer(buffer, message);

    int result = write_bytes(buffer, len, poll_id);
    if (result < 0)
    {
        return result; // error.
    }

    return 1;
}

/* Function: send_heartbeat ===============================================
 * Abstract:
 *      Send a heartbeat MAVLink messages.
 * 
 * Returns:
 *      This function will return 1 when the messages has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_heartbeat()
{
    mavlink_heartbeat_t heartbeat;
    heartbeat.type = 2; // Quadcopter
    heartbeat.autopilot = 12; // PX4
    heartbeat.base_mode = 32; // HIL mode
    heartbeat.custom_mode = 0; // None
    heartbeat.system_status = 0; // Unknown

    mavlink_message_t message;
    if(mav_sys_id < 0 || mav_comp_id < 0)
    {
        mavlink_msg_heartbeat_encode(0, 0, &message, &heartbeat);
    }
    else
    {
        mavlink_msg_heartbeat_encode(mav_sys_id, mav_comp_id, &message, &heartbeat);
    }

    return write_message(&message, get_poll_index());
}

/* Function: send_cmd_ack ===============================================
 * Abstract:
 *    Send a command acknowledge MAVLink message.
 * 
 * Parameters:
 *      cmd: the command to acknowledge.
 *      result: the acknowledgement result.
 * 
 * Returns:
 *      This function will return 1 when the messages has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_cmd_ack(uint16_t cmd, uint8_t result)
{
    mavlink_command_ack_t ack;
    ack.command = cmd;
    ack.result = result;

    mavlink_message_t message;
    mavlink_msg_command_ack_encode(mav_sys_id, mav_comp_id, &message, &ack);

    return write_message(&message, get_poll_index());
}

/* Function: send_cmd_set_mode ===============================================
 * Abstract:
 *      Sets the autopilot to HIL mode, disarmed.
 * 
 * Returns:
 *      This function will return 1 when the messages has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_cmd_set_mode()
{
    mavlink_set_mode_t cmd;
    cmd.target_system = mav_sys_id;
    cmd.base_mode = 32;

    mavlink_message_t message;
    mavlink_msg_set_mode_encode(mav_sys_id, mav_comp_id, &message, &cmd);

    return write_message(&message, get_poll_index());
}

/* Function: send_cmd_arm ===============================================
 * Abstract:
 *      Arms or disarms the autopilot.
 * 
 * Parameters:
 *      arm: the arm command.
 * 
 * Returns:
 *      This function will return 1 when the messages has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_cmd_arm(uint8_t arm)
{
    mavlink_command_long_t cmd;
    cmd.target_system = mav_sys_id;
    cmd.target_component = 0;
    cmd.command = MAV_CMD_COMPONENT_ARM_DISARM;
    cmd.confirmation = 0;
    cmd.param1 = arm;

    mavlink_message_t message;
    mavlink_msg_command_long_encode(mav_sys_id, mav_comp_id, &message, &cmd);

    return write_message(&message, get_poll_index());
}

/* Function: send_hil_messages ===============================================
 * Abstract:
 *      Send all the HIL MAVLink messages.
 * 
 * Returns:
 *      This function will return 1 when the messages has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_hil_messages(uint64_t time_usec, double q[4], double lat_lon_alt[3], double vel_e[3], double vel, double cog, double eph, double epv, int fix_type, int num_sats, double acc_b[3], double gyro[3], double mag[3], double pressure, double temperature)
{
    static uint64_t hil_state_update = 0;
    static uint64_t hil_sensor_update = 0;
    static uint64_t hil_gps_update = 0;

    int ret = 1;
    // if(hil_state_freq > 0 && (int64_t)(time_usec - hil_state_update) >= HIL_MSG_TIME_ERROR)
    // {
    //     ret = send_hil_state(time_usec, q, euler_rates, lat_lon_alt, vel_e, acc_b);
    //     hil_state_update = time_usec + (uint64_t)(1000000.0 / hil_state_freq);
    // }
    if(ret > 0 && time_usec >= time_usec_start + hil_sensor_startup_delay * 1000 && (int64_t)(time_usec - hil_sensor_update) >= HIL_MSG_TIME_ERROR)
    {
        ret = send_hil_sensors(time_usec, acc_b, gyro, mag, pressure, lat_lon_alt[2], temperature);
        hil_sensor_update = time_usec + (uint64_t)(1000000.0 / hil_sensor_freq);
    }
    if(ret > 0 && time_usec >= time_usec_start + hil_gps_startup_delay * 1000 && (int64_t)(time_usec - hil_gps_update) >= HIL_MSG_TIME_ERROR)
    {
        ret = send_hil_gps(time_usec, lat_lon_alt, vel_e, vel, cog, eph, epv, fix_type, num_sats);
        hil_gps_update = time_usec + (uint64_t)(1000000.0 / hil_gps_freq);
    }

    return ret;
}

/* Function: send_hil_state ===============================================
 * Abstract:
 *      Send the given hil_state_quaternion MAVLink message.
 * 
 * Returns:
 *      This function will return 1 when the message has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_hil_state(uint64_t time_usec, double q[4], double euler_rates[3], double lat_lon_alt[3], double vel_e[3], double acc_b[3])
{
    mavlink_hil_state_quaternion_t hil_state;
    hil_state.time_usec = time_usec;
    hil_state.attitude_quaternion[0] = (float)q[0];
    hil_state.attitude_quaternion[1] = (float)q[1];
    hil_state.attitude_quaternion[2] = (float)q[2];
    hil_state.attitude_quaternion[3] = (float)q[3];
    hil_state.rollspeed = (float)euler_rates[0];
    hil_state.pitchspeed = (float)euler_rates[1];
    hil_state.yawspeed = (float)euler_rates[2];
    hil_state.lat = (int32_t)(lat_lon_alt[0] * 1e7);
    hil_state.lon = (int32_t)(lat_lon_alt[1] * 1e7);
    hil_state.alt = (int32_t)(lat_lon_alt[2] * 1e3);
    hil_state.vx = (int16_t)(vel_e[0]*100);
    hil_state.vy = (int16_t)(vel_e[1]*100);
    hil_state.vz = (int16_t)(vel_e[2]*100);
    hil_state.xacc = (int16_t)(acc_b[0] * 1000.0f / 9.81f);
    hil_state.yacc = (int16_t)(acc_b[1] * 1000.0f / 9.81f);
    hil_state.zacc = (int16_t)(acc_b[2] * 1000.0f / 9.81f);
    hil_state.true_airspeed = (uint16_t)(sqrt(pow(vel_e[0], 2) + pow(vel_e[1], 2) + pow(vel_e[2], 2)) * 100);

    mavlink_message_t message;
    mavlink_msg_hil_state_quaternion_encode(mav_sys_id, mav_comp_id, &message, &hil_state);

    return write_message(&message, get_poll_index());
}

/* Function: send_hil_sensors ===============================================
 * Abstract:
 *      Send the given hil_sensor MAVLink message.
 * 
 * Returns:
 *      This function will return 1 when the message has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_hil_sensors(uint64_t time_usec, double acc_b[3], double gyro[3], double mag[3], double pressure, double alt, double temperature)
{
    mavlink_hil_sensor_t hil_sensor;
    hil_sensor.time_usec = time_usec;
    hil_sensor.xacc = (float)(acc_b[0]);
    hil_sensor.yacc = (float)(acc_b[1]);
    hil_sensor.zacc = (float)(acc_b[2]);
    hil_sensor.xgyro = (float)(gyro[0]);
    hil_sensor.ygyro = (float)(gyro[1]);
    hil_sensor.zgyro = (float)(gyro[2]);
    hil_sensor.xmag = (float)(mag[0]);
    hil_sensor.ymag = (float)(mag[1]);
    hil_sensor.zmag = (float)(mag[2]);
    hil_sensor.abs_pressure = (float)(pressure * 0.01);
    hil_sensor.pressure_alt = (float)(alt * 1e3);
    hil_sensor.temperature = (float)temperature;
    hil_sensor.fields_updated = 0x1FFF;

    mavlink_message_t message;
    mavlink_msg_hil_sensor_encode(mav_sys_id, mav_comp_id, &message, &hil_sensor);

    return write_message(&message, get_poll_index());
}

/* Function: send_hil_gps ===============================================
 * Abstract:
 *      Send the given hil_gps MAVLink message.
 * 
 * Returns:
 *      This function will return 1 when the message has been sent.
 *      If an error ocurred, -1 will be returned.
 */
int send_hil_gps(uint64_t time_usec, double lat_lon_alt[3], double vel_e[3], double vel, double cog, double eph, double epv, int fix_type, int num_sats)
{
    mavlink_hil_gps_t hil_gps;
    hil_gps.time_usec = time_usec;
    hil_gps.fix_type = (uint8_t)fix_type;
    hil_gps.lat = (int32_t)(lat_lon_alt[0] * 1e7);
    hil_gps.lon = (int32_t)(lat_lon_alt[1] * 1e7);
    hil_gps.alt = (int32_t)(lat_lon_alt[2] * 1e3);
    hil_gps.eph = (uint16_t)(eph * 100);
    hil_gps.epv = (uint16_t)(epv * 100);
    hil_gps.vn = (int16_t)(vel_e[0] * 100);
    hil_gps.ve = (int16_t)(vel_e[1] * 100);
    hil_gps.vd = (int16_t)(vel_e[2] * 100);
    hil_gps.vel = (uint16_t)(vel * 100);
    hil_gps.cog = (uint16_t)(cog * 100);
    hil_gps.satellites_visible = (uint8_t)num_sats;

    mavlink_message_t message;
    mavlink_msg_hil_gps_encode(mav_sys_id, mav_comp_id, &message, &hil_gps);

    return write_message(&message, get_poll_index());
}

void get_thrust_commands_force(double thrust_commands[4])
{
    thrust_commands[0] = (max_thrust / 4.0) * ((!hover_reached) ? thrust_hover : (double)hil_actuators.controls[0]);//2
    thrust_commands[1] = (max_thrust / 4.0) * ((!hover_reached) ? thrust_hover : (double)hil_actuators.controls[1]);//0
    thrust_commands[2] = (max_thrust / 4.0) * ((!hover_reached) ? thrust_hover : (double)hil_actuators.controls[2]);//3
    thrust_commands[3] = (max_thrust / 4.0) * ((!hover_reached) ? thrust_hover : (double)hil_actuators.controls[3]);//1
}