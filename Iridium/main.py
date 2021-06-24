""" ############################################################################
#    Main PROGRAM OF FZ, Written by T.K.             Last updated : Jun. 7, 2021
#
# -Coding Rules
#   - fz.prt function is only used here except for the debugging/error purpose
#
# -Notes
#   The parameters specified in param.txt are as follows
#   - T    : Sampling interval (second)
#   - dt   : Sampling interval (second)
#   - Nw   : The window size of the segment for the spectral analysis.
#   - Ncut : The index of Cut-off Frequency for the frequency spectral Pzz.
#   - imu_i : imu_i switches the IMU for logging as follows
#       - 1 : BNO055
#       - 2 : BNO08x
#       - 3 : LSM6DS33
#       - 4 : ISM330DXCX
#   - gps_i : GPS data (GPSRMC & GPSGGA) are logged before logging IMU.
#   - pix_val : The intensity of Neopixel (0-255)
#
# -Add iridium code by T.Katsuno : Jun. 24, 2021
#   ############################################################################
"""
import gc
import fz
import time as t

# for RockBlock iridium
import struct
import board
import adafruit_rockblock

# Initizalize
#   - 1. Import parameters
T, dt, Nw, Ncut, frq, imu_i, gps_i, psd_i, pix_val = fz.ior_prm()
#   - 2. Setup each sensors (IMU, GPS, RTC, SD, Neopixel)
set_rtc = True
set_t = (2021, 6, 24, 17, 00, 0, 0, -1, -1)
pixels, i2c, imu, gps, Dir_Out = fz.dvc_ini(imu_i, gps_i, set_rtc, set_t)

# RockBlock setup
uart = board.UART()
uart.baudrate = 19200
rb = adafruit_rockblock.RockBlock(uart)

while True:

    # Start measurement every calendar 15min
    while not (fz.rtc_now(i2c).tm_min % 15 == 0) :
        fz.npx_blk2(1, pixels)

    # Initizalize for each measurement
    time_now, ts, Fil_Log = fz.log_ini(i2c, pixels, pix_val)
    fz.prt(Fil_Log, "\n****************************************")
    fz.prt(Fil_Log, "\n\nPRORGRAM FZ start ")
    fz.prt(Fil_Log, "\nT={:} (s)\t dt={:} (s)\t Nw={:} ".format(T, dt, Nw))
    fz.prt(Fil_Log, "\nLog file is {}".format(Fil_Log))
    fz.prt(Fil_Log, "\n\t Measurement started on {} (yyyymmdd_HHMM)".format(time_now))

    # Logging GPS (only if gps_i == 1)
    if gps_i == 1 :
        fz.prt(Fil_Log, "\nRecording GPS data ... ")
        while not gps.has_fix:
            gps.update()
        for i in range(2) :
            fz.prt(Fil_Log, "\n\t {}".format(gps.readline()))

    # Logging IMU (the sensor can be switched by imu_i )
    fz.prt(Fil_Log, "\nRecording IMU data ... {:6.1f} sec".format(t.monotonic()-ts))
    Fil_raw, time_now, N_l = fz.iow_imu(Dir_Out, time_now, i2c, imu, imu_i, pixels)
    fz.prt(Fil_Log, "\n\t Logging finished. {1:5d} lines -> {0:}".format(Fil_raw, N_l))

    # Calculate wave statistics from IMU data
    if psd_i == 1 :
        fz.prt(Fil_Log, "\n Data Analysis ... {:6.1f} sec".format(t.monotonic()-ts))
        fz.prt(Fil_Log, "\n\tLow pass filtering : {:6.1f} sec".format(t.monotonic()-ts))
        Fil_low, N_low = fz.flt_low(Fil_raw, N_l, dt)

        fz.prt(Fil_Log, "\n\tEstimates psd　    : {:6.1f} sec".format(t.monotonic()-ts))
        Pxx, Pyy, Pzz, Qxz, Qyz, Cxy = fz.psd_seg(Fil_low, N_low)

        fz.prt(Fil_Log, "\n\tCal wave stats.　  : {:6.1f} sec".format(t.monotonic()-ts))
        Hs, Hsx, Hsy, Fp, Fpx, Fpy, Mwd, Mds = fz.blk_sta(Pxx, Pyy, Pzz, Qxz, Qyz, Cxy)

        fz.prt(Fil_Log, "\n\t\tHs={:5.3f}(m), Fp={:4.2f}(Hz)".format(Hs, Fp))
        fz.prt(Fil_Log, "\n\t\tHs_x={:5.3f}(m), Fp_x={:4.2f}(Hz)".format(Hsx, Fpx))
        fz.prt(Fil_Log, "\n\t\tHs_y={:5.3f}(m), Fp_z={:4.2f}(Hz)".format(Hsy, Fpy))

        # Rockblock iridium transfer
        data = struct.pack("f", Hs)
        data += struct.pack("f", Fp)
        if gps_i == 1 :
            while not gps.has_fix :
                gps.update()
            for i in range(2) :
                data += struct.pack("f", gps.latitude)
                data += struct.pack("f", gps.longitude)
        retry = 1
        rb.data_out = data
        status = rb.satellite_transfer()
        print(retry, status)
        t.sleep(10)
        while status[0] > 8:
            retry += 1
            rb.data_out = data
            status = rb.satellite_transfer()
            print(retry, status)
            t.sleep(10)

    # Finished ... wait for the next loop that starts just every T seconds
    #   - 1. the memory is checked.
    #   - 2. Neopixel turns to be white.
    fz.prt(Fil_Log, "\ngc.mem_free() = {:6.1f}".format(gc.mem_free()))
    fz.prt(Fil_Log, "\n******************")
    pixels[0] = (pix_val, pix_val, pix_val)
