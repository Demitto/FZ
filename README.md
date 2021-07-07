# FZ                                         Last updated : Jul. 7, 2021

This is repository for circuitpython codes for using IMUs produced by adafruit.

There are the following three main files 
  1. main.py  
  2. param.txt
  3. ./lib/fz.py

  - The main.py is designed to overview the whole processes, while each process is written in fz.py. The runtime parameters can be specified by param.txt. Description for each parameter can be found in the head of main.py and also at the end of this README.md.

NeoPixel behavior 
  1. blink quickly (blue/green) once the board is powered.
  2. IMU calibration status from red-orange-white-blue. 
  3. GPS: Red is for initializing / wake up.  Orange is for searching fix.
  4. blink every .5 seconds (blue/green) before the logging start every T seconds
  5. blue without blinking when the data are being logged. 

  - Note: One needs to calibrate before the logging starts if cal_i is set to 1.

Getting started with adafruit-feather with circuitpy 
  1. Double click the reset button following the official description  https://learn.adafruit.com/adafruit-feather-sense/circuitpython-on-feather-sense
  2. Copy bno055.zip from this git-hub repository for the jump-start. 
  3. Play with the feather & IMUs. 


****************************************************************************
*    Main PROGRAM OF FZ, Written by T.K.             Last updated : Jul. 5, 2021
*
* -Coding Rules
*   - fz.prt function is only used here except for the debugging/error purpose
*
* -Notes
*   The parameters specified in param.txt are as follows
*   - T       : Measurement Cycle [s]
*   - T_imu   : IMU logging time [s] (Note: T_imu < T)
*   - Hz1     : IMU Sampling Frequency 1 [Hz]
*   - Hz2     : IMU Output interval 2 [Hz]
*   - Nw      : The window size of the segment for the spectral analysis.
*   - Ncut    : The index of Cut-off Frequency for the frequency spectral Pzz.
*   - imu_i   : if 1, IMU data are logged.
*   - cal_i   : if 1, IMU is calibrated before the first measurement (imu_i should be 1)
*   - psd_i   : if 1, PSD and bulk wave statistics are calculated (imu_i should be 1)
*   - gps_i   : if 1, GPS data aree logged
*   - sen_i   : if 1, Feather Sense data (Air Temp., Air Pres., Noise Level) are logged
*   - pix_val (pval) : The intensity of Neopixel (0-255)
****************************************************************************
