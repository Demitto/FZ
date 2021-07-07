# FZ                                         Last updated : Jul. 7, 2021

This is repository for circuitpython codes for using IMUs produced by adafruit.

There are the following three main files 
  1. main.py  
  2. param.txt
  3. ./lib/fz.py

  - The main.py is designed to overview the whole processes, while each process is written in fz.py. The runtime parameters can be specified by param.txt. Description for each parameter can be found in the head of main.py.

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
