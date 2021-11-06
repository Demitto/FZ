# FZ                                         Last updated : Nov. 6, 2021

This is repository for circuitpython codes for using IMUs produced by adafruit.

There are the following three main files 
  1. main.py  
  2. param.txt
  3. ./lib/fz.py
  - The main.py is designed to overview the whole processes, while each process is written in fz.py. The runtime parameters can be specified by param.txt. Description for each parameter can be found in the head of main.py.

NeoPixel behavior 
  1. blink quickly (blue/green) once the board is powered.
  2. blue is used for IMU: blinking until IMU logging start.
  3. green is used for GPS: blinking until GPS fix is obtained.
  - Note: Red-Blue blinking is used for IMU calibration (only happen with cal_i = 1)
