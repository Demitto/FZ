""" ############################################################################
#    Main PROGRAM OF FZ, Written by T.K.             Last updated : Jun. 29, 2021
#
# -Coding Rules
#   - fz.prt function is only used here except for the debugging/error purpose
#
# -Notes
#   The parameters specified in param.txt are as follows
#   - T    : Sampling interval (second)
#   - 1/Hz1   : Sampling interval 1 (second)
#   - 1/Hz2   : Sampling interval 2 (second)
#   - Nw   : The window size of the segment for the spectral analysis.
#   - Ncut : The index of Cut-off Frequency for the frequency spectral Pzz.
#   - imu_i : imu_i switches the IMU for logging as follows
#       - 1 : BNO055
#       - 2 : BNO08x
#       - 3 : LSM6DS33
#       - 4 : ISM330DXCX
#   - gps_i : GPS data (GPSRMC & GPSGGA) are logged before logging IMU.
#   - pix_val (pval) : The intensity of Neopixel (0-255)
#   ############################################################################
"""
import gc
from math import sqrt as sqrt
import fz
import time as t
import struct

# Initizalize
#   - 0. Optional: Set RTC if set_rtc == True
sett = False
rt = (2021, 6, 30, 20, 41, 0, 0, -1, -1)
#   - 1. Import parameters
T, Hz1, Hz2, Nw, Ncut, frq, imu_i, cal_i, gps_i, sen_i, psd_i, iri_i, pval = fz.ior_prm()
#   - 2. Setup each sensors (IMU, GPS, RTC, SD, Neopixel)
pix, i2c, rtc, imu, gps, bmp, mic, micv, hmd, rb = \
    fz.dvc_ini(imu_i, cal_i, gps_i, sen_i, iri_i, pval, sett, rt)

# Measurement Loop
while True:

    # Start measurement every calendar int(T/60) minutes
    if True :
        while not (rtc.datetime.tm_min % int(T/60) == 0) :
            fz.npx_blk(1, pix, pval)

    # Initizalize for each measurement
    time_now, ts, Fil_Log = fz.log_ini(i2c, pix, pval)
    fz.prt(Fil_Log, "\n****************************************")
    fz.prt(Fil_Log, "\n\nPRORGRAM FZ start ")
    fz.prt(Fil_Log, "\nT={:} (s)\t dt={:} (s) \t Nw={:} ".format(T, 1/Hz1, Nw))
    fz.prt(Fil_Log, "\nLog file is {}".format(Fil_Log))
    fz.prt(Fil_Log, "\n\t Measurement started on {} (yyyymmdd_HHMM)".format(time_now))

    # Logging GPS (only if gps_i == 1)
    if gps_i == 1 :
        fz.prt(Fil_Log, "\nRecording GPS data ... ")
        while not gps.has_fix:
            gps.update()
        for i in range(2) :
            fz.prt(Fil_Log, "\n\t {}".format(str(gps.readline(), "ascii").strip()))
        if(gps.timestamp_utc.tm_year > 2000) :
            set_time = (gps.timestamp_utc.tm_year, gps.timestamp_utc.tm_mon,
                        gps.timestamp_utc.tm_mday, gps.timestamp_utc.tm_hour,
                        gps.timestamp_utc.tm_min, gps.timestamp_utc.tm_sec, 0, -1, -1)
            rtc.datetime = t.struct_time(set_time)
        t_n = rtc.datetime
        now = "{:04d}{:02d}{:02d}_{:02d}{:02d}{:02d}".format(
            t_n.tm_year, t_n.tm_mon, t_n.tm_mday, t_n.tm_hour, t_n.tm_min, t_n.tm_sec)
        fz.prt(Fil_Log, "\n\t The RTC time is set to {} (yyyymmdd_HHMMSS)".format(now))

    # Logging BMP280 temperature and barometric pressure(only if sen_i == 1)
    if sen_i == 1 :
        fz.prt(Fil_Log, "\nRecording Feather Sense data ... ")
        fz.prt(Fil_Log, "\n\tAir Temperaure : {:.2f} C".format(bmp.temperature))
        fz.prt(Fil_Log, "\n\tAir Pressure   : {:.2f} hPa".format(bmp.pressure))
        fz.prt(Fil_Log, "\n\tAir Humidity   : {:.2f} %".format(hmd.relative_humidity))
        minb = int(sum(micv) / len(micv))
        fz.prt(Fil_Log, "\n\tSound level    : {} ".format(
               int(sqrt(sum(float(x - minb) * (x - minb) for x in micv) / len(micv)))))

    # Logging IMU (the sensor can be switched by imu_i )
    fz.prt(Fil_Log, "\nRecording IMU data ... {:6.1f} sec".format(t.monotonic()-ts))
    Fil_low, time_now, N_l = fz.iow_imu(time_now, i2c, imu, imu_i, pix)
    fz.prt(Fil_Log, "\n\tLogging finished. {1:5d} lines -> {0:}".format(Fil_low, N_l))

    # Calculate wave statistics from IMU data
    if psd_i == 1 :
        fz.prt(Fil_Log, "\n Data Analysis ... {:6.1f} sec".format(t.monotonic()-ts))
        fz.prt(Fil_Log, "\n\tEstimates psd　    : {:6.1f} sec".format(t.monotonic()-ts))
        gc.collect()
        fz.prt(Fil_Log, "\n\tgc.mem_free() = {:6.1f}".format(gc.mem_free()))
        Pxx, Pyy, Pzz, Qxz, Qyz, Cxy = fz.psd_seg(Fil_low, N_l)
        fz.prt(Fil_Log, "\n\tCal wave stats.　  : {:6.1f} sec".format(t.monotonic()-ts))
        Hs, Hsx, Hsy, Fp, Fpx, Fpy, Mwd, Mds = fz.blk_sta(Pxx, Pyy, Pzz, Qxz, Qyz, Cxy)
        fz.prt(Fil_Log, "\n\t\tHs={:5.3f}(m), Fp={:4.2f}(Hz)".format(Hs, Fp))
        fz.prt(Fil_Log, "\n\t\tHs_x={:5.3f}(m), Fp_x={:4.2f}(Hz)".format(Hsx, Fpx))
        fz.prt(Fil_Log, "\n\t\tHs_y={:5.3f}(m), Fp_z={:4.2f}(Hz)".format(Hsy, Fpy))

    #Transfer data using iridium
    if iri_i == 1 :
        fz.prt(Fil_Log, "\n Iridium Transfer ... {:6.1f} sec".format(t.monotonic()-ts))
        data = struct.pack("f", Hs)
        data += struct.pack("f", Fp)
        if gps_i == 1 :
            while not gps.has_fix :
                gps.update()
            data += struct.pack("f", gps.latitude)
            data += struct.pack("f", gps.longitude)
        fz.iri(rb, data)
        fz.prt(Fil_Log, "\n\t data={}".format(data))

    # Finished ... wait for the next loop that starts just every T seconds
    #   - 1. the memory is checked.
    #   - 2. Neopixel turns to be white.
    fz.prt(Fil_Log, "\ngc.mem_free() = {:6.1f} bef. gc.collect()".format(gc.mem_free()))
    gc.collect()
    fz.prt(Fil_Log, "\ngc.mem_free() = {:6.1f} aft. gc.collect()".format(gc.mem_free()))
    fz.prt(Fil_Log, "\n******************")
    pix[0] = (pval, pval, pval)
