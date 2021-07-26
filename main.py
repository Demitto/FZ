""" ############################################################################
#    Main PROGRAM OF FZ, Written by T.K.             Last updated : Jul. 5, 2021
#
# -Coding Rules
#   - fz.prt function is only used here except for the debugging/error purpose
#
# -Notes
#   The parameters specified in param.txt are as follows
#   - T       : Measurement Cycle [s]
#   - T_imu   : IMU logging time [s] (Note: T_imu < T)
#   - Hz1     : IMU Sampling Frequency 1 [Hz]
#   - Hz2     : IMU Output interval 2 [Hz]
#   - Nw      : The window size of the segment for the spectral analysis.
#   - Ncut    : The index of Cut-off Frequency for the frequency spectral Pzz.
#   - imu_i   : if 1, IMU data are logged.
#   - cal_i   : if 1, IMU is calibrated before the first measurement (imu_i should be 1)
#   - psd_i   : if 1, PSD and bulk wave statistics are calculated (imu_i should be 1)
#   - gps_i   : if 1, GPS data aree logged
#   - sen_i   : if 1, Feather Sense data (Air Temp., Air Pres., Noise Level) are logged
#   - pix_val (pval) : The intensity of Neopixel (0-255)
#   ############################################################################
"""
import gc
from math import sqrt as sqrt
import fz
import time as t

# 1. Initizalize
#   - 0. Optional: Set RTC if set_rtc == True
sett = False
rt = (2021, 7, 2, 13, 15, 5, 0, -1, -1)
#    1-1. Import parameters
T, T_imu, Hz1, Hz2, Nw, Ncut, frq, imu_i, cal_i, gps_i, sen_i, psd_i, pval = \
    fz.ior_prm()
#    1-2. Setup each sensors (IMU, GPS, RTC, SD, Neopixel)
pix, i2c, rtc, imu, gps, bmp, mic, micv, hmd = \
    fz.dvc_ini(T, imu_i, cal_i, gps_i, sen_i, pval, sett, rt)

# 2. Measurement Loop
while True:

    # 2-0. Start measurement every calendar int(T/60) minutes
    if True :
        fz.imu_chg(imu, 0)
        while not (rtc.datetime.tm_min % int(T/60) == 0) :
            fz.npx_blk(1, pix, [0, 0, 0], [0, 0, pval])
        fz.imu_chg(imu, 1)
    # Green !
    pix[0] = (0, pval, 0)

    # 2-1. Initizalize for each measurement
    time_now, ts, Fil_Log = fz.log_ini(i2c, rtc, pix, pval)
    fz.prt(Fil_Log, "\n************************************************************")
    fz.prt(Fil_Log, "\n\t\tFZ Measurement Log")
    fz.prt(Fil_Log, "\n************************************************************")
    fz.prt(Fil_Log, "\nMeasurement starts ... \t {:6.1f} [s]".format(t.monotonic()-ts))
    fz.prt(Fil_Log, "\n\tThe start date_time is {} (yyyymmdd_HHMM)".format(time_now))
    fz.prt(Fil_Log, "\n\tLog file is {}".format(Fil_Log))
    fz.prt(Fil_Log, "\n\tIMU data file is IMU_BNO055_{0}.txt".format(time_now))
    # fz.prt(Fil_Log, "\n\tInput parameters are as follows \r\n")
    # with open("param.txt", "r") as fp:
    #    for line in fp:
    #        fz.prt(Fil_Log, "\t\t{}".format(line))

    # 2-2.Logging GPS (only if gps_i == 1)
    if gps_i == 1 :
        fz.prt(Fil_Log, "\nRecording GPS data ...\t{:6.1f}[s]".format(t.monotonic()-ts))
        fz.gps_ini(i2c, gps, rtc, pix, pval)
        if((t.monotonic()-ts) > (T - T_imu)*1.5) :
            fz.prt(Fil_Log, "\nRecording GPS failed ...")
            t.sleep(max(T - (t.monotonic()-ts) - 10 , 0))
            continue
        nmea1, nmea2, lat, lon = fz.gps_log(gps, T - (t.monotonic()-ts) - 10, pix, pval)
        fz.prt(Fil_Log, "\n\tlatitude = {:9.4f}".format(lat))
        fz.prt(Fil_Log, "\n\tlongitude = {:9.4f}".format(lon))
        fz.prt(Fil_Log, "\n\tNMEA output (RMC and GGA)")
        fz.prt(Fil_Log, "\n\t {}\n".format(nmea1))
        fz.prt(Fil_Log, "\n\t {}\n".format(nmea2))

    # 2-3. Logging BMP280 temperature and barometric pressure(only if sen_i == 1)
    if sen_i == 1 :
        fz.prt(Fil_Log, "\nRecording Air data ...\t{:6.1f}[s]".format(t.monotonic()-ts))
        fz.prt(Fil_Log, "\n\tAir Temperaure : {:.2f} C".format(bmp.temperature))
        fz.prt(Fil_Log, "\n\tAir Pressure   : {:.2f} hPa".format(bmp.pressure))
        fz.prt(Fil_Log, "\n\tAir Humidity   : {:.2f} %".format(hmd.relative_humidity))
        minb = int(sum(micv) / len(micv))
        fz.prt(Fil_Log, "\n\tSound level    : {} ".format(
               int(sqrt(sum(float(x - minb) * (x - minb) for x in micv) / len(micv)))))

    # 2-4. Logging IMU (the sensor can be switched by imu_i )
    fz.prt(Fil_Log, "\n\nRecording IMU data ...\t{:6.1f} [s]".format(t.monotonic()-ts))
    # Blue !
    pix[0] = (0, 0, pval)
    Fil_low, time_now, N_l = \
        fz.iow_imu(time_now, T, T_imu, Hz1, Hz2, i2c, imu, imu_i, pix)
    fz.prt(Fil_Log, "\n\tLogging finished ...")
    fz.prt(Fil_Log, "\n\t\t The number of lines : {:5d}".format(N_l))
    fz.prt(Fil_Log, "\n\t\t ---> {}\n".format(Fil_low))

    # 3. Calculate wave statistics from IMU data
    if psd_i == 1 :
        fz.prt(Fil_Log, "\nIMU Data Analysis ...\t{:6.1f}[s]".format(t.monotonic()-ts))
        fz.prt(Fil_Log, "\n\tEstimates PSD ...")
        gc.collect()
        # fz.prt(Fil_Log, "\n\tgc.mem_free() = {:6.1f}".format(gc.mem_free()))
        Pxx, Pyy, Pzz, Qxz, Qyz, Cxy = fz.psd_seg(Fil_low, N_l, Nw, frq)
        fz.prt(Fil_Log, "\n\tCalculate bulk wave statistics ...")
        Hs, Hsx, Hsy, Fp, Fpx, Fpy, Mwd, Mds = \
            fz.blk_sta(Ncut, frq, Pxx, Pyy, Pzz, Qxz, Qyz, Cxy)
        fz.prt(Fil_Log, "\n\t\tHs={:5.3f}(m), Fp={:4.2f}(Hz)".format(Hs, Fp))
        fz.prt(Fil_Log, "\n\t\tHs_x={:5.3f}(m), Fp_x={:4.2f}(Hz)".format(Hsx, Fpx))
        fz.prt(Fil_Log, "\n\t\tHs_y={:5.3f}(m), Fp_z={:4.2f}(Hz)".format(Hsy, Fpy))

    # 4. Finished
    gc.collect()
    fz.prt(Fil_Log, "\n\nMeasurement finished :\t{:6.1f} [s]".format(t.monotonic()-ts))
    fz.prt(Fil_Log, "\n************************************************************")
