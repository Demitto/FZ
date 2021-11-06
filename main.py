""" ############################################################################
#    Main PROGRAM OF FZ, Written by T.K.             Last updated : Nov. 2, 2021
#
# -Coding Rules
#   - The output to the serial console is only made in this main.py. More specifically,
#     fz.prt function is only used here except for the debugging/error purpose
#
# -Notes
#   The parameters specified in param.txt are as follows
#   - T1      : Measurement Cycle [s]
#   - T2      : IMU logging time [s] (Note: T_imu < T)
#   - Hz1     : IMU Sampling Frequency 1 [Hz]
#   - Hz2     : IMU Output interval 2 [Hz]
#   - Nw      : The window size of the segment for the spectral analysis.
#   - Nc      : The index of Cut-off Frequency for the frequency spectral Pzz.
#   - imu_i   : if 1, IMU data are logged.
#   - cal_i   : if 1, IMU is calibrated before the first measurement (imu_i should be 1)
#   - psd_i   : if 1, PSD and bulk wave statistics are calculated (imu_i should be 1)
#   - gps_i   : if 1, GPS data aree logged
#   - sen_i   : if 1, Feather Sense data (Air Temp., Air Pres., Noise Level) are logged
#   - iri_i   : if 1, Data is transmitted via Iridium Satellite Burst Data (Iridium-SBD)
#   ############################################################################
"""
# 0. Module import
if True :
    import gc
    from math import sqrt as sqrt
    import time as t
    import fz
    import supervisor

# 1. Initizalization
if True :
    #   - 0. Optional: Set RTC if set_rtc == True
    sett = False
    rt = (2021, 7, 2, 13, 15, 5, 0, -1, -1)
    #    1-1. Import parameters
    T1, T2, Hz1, Hz2, Nw, Nc, frq, imu_i, cal_i, gps_i, sen_i, psd_i, iri_i, wdt_i = fz.prm()
    #    1-2. Setup each sensors (IMU, GPS, RTC, SD, Neopixel)
    npx, i2c, rtc, imu, mgn, gps, bmp, mic, micv, hmd, gps_i, rb = \
        fz.ini_all(imu_i, cal_i, gps_i, sen_i, iri_i, sett, rt)
    #    1-3. Setup watch dog timer
    if wdt_i == 1 :
        wdt = fz.ini_wdt(T1*2)

# 2. Measurement Loop
try :

    while True:
        # 2-1. Initizalize for each measurement
        time_now, ts, Fil_Log = fz.ini_log(T1, imu, rtc, npx)
        fz.prt(Fil_Log, "\n***********************************************************")
        fz.prt(Fil_Log, "\n\t\tFZ Measurement Log")
        fz.prt(Fil_Log, "\n***********************************************************")
        fz.prt(Fil_Log, "\n Measurement started on {} (yyyymmdd_HHMM)".format(time_now))
        fz.prt(Fil_Log, "\n\tLog file is {}".format(Fil_Log))
        fz.prt(Fil_Log, "\n\tIMU data file is IMU_BNO055_{0}.txt".format(time_now))
        # fz.prt(Fil_Log, "\n\tInput parameters are as follows \r\n")
        # with open("param.txt", "r") as fp:
        #    for line in fp:
        #        fz.prt(Fil_Log, "\t\t{}".format(line))
        # 2-2.Logging GPS (only if gps_i == 1)
        if gps_i == 1 :
            fz.prt(Fil_Log, "\nGPS data ...\t{:6.1f}[s]".format(t.monotonic()-ts))
            nmea1, nmea2, lat, lon = fz.log_gps(gps, rtc, npx)
            fz.prt(Fil_Log, "\n\tlatitude = {:9.4f}".format(lat))
            fz.prt(Fil_Log, "\n\tlongitude = {:9.4f}".format(lon))
            fz.prt(Fil_Log, "\n\tNMEA output (RMC and GGA)")
            fz.prt(Fil_Log, "\n\t {}\n".format(nmea1))
            fz.prt(Fil_Log, "\n\t {}\n".format(nmea2))
        # 2-3. Logging BMP280 temperature and barometric pressure(only if sen_i == 1)
        if sen_i == 1 :
            fz.prt(Fil_Log, "\nAir data ...\t {:6.1f}[s]".format(t.monotonic()-ts))
            fz.prt(Fil_Log, "\n\tAir Temp.  : {:.2f} C".format(bmp.temperature))
            fz.prt(Fil_Log, "\n\tAir Pres.  : {:.2f} hPa".format(bmp.pressure))
            fz.prt(Fil_Log, "\n\tAir Humi.  : {:.2f} %".format(hmd.relative_humidity))
            minb = int(sum(micv) / len(micv))
            fz.prt(Fil_Log, "\n\tSound level    : {} ".format(
                int(sqrt(sum(float(x - minb) * (x - minb) for x in micv) / len(micv)))))
        # 2-4. Logging IMU (the sensor can be switched by imu_i )
        if imu_i == 1 :
            fz.prt(Fil_Log, "\n\nLogging IMU ...\t{:6.1f} [s]".format(t.monotonic()-ts))
            Fil_low, Nl = fz.log_imu(time_now, T2, Hz1, Hz2, imu, mgn, npx)
            fz.prt(Fil_Log, "\n\tLogging finished ...")
            fz.prt(Fil_Log, "\n\t\t The number of lines : {:5d}".format(Nl))
            fz.prt(Fil_Log, "\n\t\t ---> {}\n".format(Fil_low))
        # 3. Calculate wave statistics from IMU data
        if psd_i == 1 :
            fz.prt(Fil_Log, "\nSpec. Analysis ...\t{:6.1f}[s]".format(t.monotonic()-ts))
            fz.prt(Fil_Log, "\n\tEstimates PSD ...")
            gc.collect()
            # fz.prt(Fil_Log, "\n\tgc.mem_free() = {:6.1f}".format(gc.mem_free()))
            Pxx, Pyy, Pzz, Qxz, Qyz, Cxy = fz.cal_spc(Fil_low, Nl, Nw, frq)
            fz.prt(Fil_Log, "\n\tCalculate bulk wave statistics ...")
            Hs, Hsx, Hsy, Fp, Fpx, Fpy, Mwd, Mds = \
                fz.cal_blk(Nc, frq, Pxx, Pyy, Pzz, Qxz, Qyz, Cxy)
            fz.prt(Fil_Log, "\n\t\tHs={:5.3f}(m), Fp={:4.2f}(Hz)".format(Hs, Fp))
            fz.prt(Fil_Log, "\n\t\tHs_x={:5.3f}(m), Fp_x={:4.2f}(Hz)".format(Hsx, Fpx))
            fz.prt(Fil_Log, "\n\t\tHs_y={:5.3f}(m), Fp_z={:4.2f}(Hz)".format(Hsy, Fpy))

        # 4. Transfer data using iridium(only if iri_i == 1)
        if iri_i == 1 :
            fz.prt(Fil_Log, "\n\nIridium SBD...\t{:6.1f} [s]".format(t.monotonic()-ts))
            rb_sent, Nr = fz.iri(rb, Hs, Fp, Pzz, lat, lon, npx, Fil_Log)
            fz.prt(Fil_Log, "\n\tTrasfrer data =\t{}".format(rb_sent))
            fz.prt(Fil_Log, "\n\tNumber of Retry = \t{}".format(Nr))

        # 5. Finished
        gc.collect()
        fz.prt(Fil_Log, "\n\nFinished after {:6.1f} [s]".format(t.monotonic()-ts))
        fz.prt(Fil_Log, "\n**********************************************************")
        wdt.feed()
# Reboot when error
except OSError :
    supervisor.reload()
