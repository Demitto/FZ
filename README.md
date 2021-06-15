# Feather-MZ                                         Last updated : Jun. 7, 2021

This is repository for circuitpy codes for using IMUs prduced by adafruit.

There are the following three main files 
  1. main.py  
  2. param.txt
  3. ./lib/fz.py
The main.py is designed to overview the whole processes, while each process is written in fz.py. The runtime parameters can be specified by param.txt. Each parameter can be described in the head of main.py and at the end of this README.md.

NeoPixel behavior 
  1. blink quickly (blue/green) once the board is powered.
  2. calibration status from red-orange-whie-blue. The 
  3. blink every .5 seconds (blue) before the logging start every 15min (??:00, ??:15, ??:30, ??:45)
  4. blue without blinking when the data are being logged. 
Note: One needs to calibrate before the logging starts. 

############################################################################

Main PROGRAM OF FZ, Written by T.K.             Last updated : Jun. 7, 2021

 -Coding Rules
   - fz.prt function is only used here except for the debugging/error purpose

 -Notes
   The parameters specified in param.txt are as follows
   - T    : Sampling interval (second)
   - dt   : Sampling interval (second)
   - Nw   : The window size of the segment for the spectral analysis.
   - Ncut : The index of Cut-off Frequency for the frequency spectral Pzz.
   - imu_i : imu_i switches the IMU for logging as follows
       - 1 : BNO055
       - 2 : BNO08x
       - 3 : LSM6DS33
       - 4 : ISM330DXCX
   - gps_i : if 1 : GPS data (GPSRMC & GPSGGA) are logged before logging IMU.
   - psd_i : if 1 : Spectral Analysis is conducted.
   - pix_val : The intensity of Neopixel (0-255)
   
############################################################################
