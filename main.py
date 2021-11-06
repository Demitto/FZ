""" ############################################################################
#    Main PROGRAM OF FZ, Written by T.K.             Last updated : Nov. 5, 2021
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
# 0. Initizalization
# import gc
import time as t
import fz
import storage
import alarm
from microcontroller import watchdog as wdt
from watchdog import WatchDogMode as wdm
#   -0-1. Optional: Set RTC if set_rtc == True
set_rtc = False
set_time = (2021, 7, 2, 13, 15, 5, 0, -1, -1)
#    0-2. Import parameters
T1, T2, Hz1, Hz2, Nw, Nc, frq, imu_i, cal_i, gps_i, sen_i, psd_i, iri_i = fz.prm()
#    0-3. Setup
npx, i2c, rtc = fz.ini_sp_(cal_i, set_rtc, set_time)
#    0-4. Watchdog
wdt.timeout = T1 * 2
wdt.mode = wdm.RESET
wdt.feed()

try :
    # 1. Waiting
    while not (rtc.datetime.tm_min % max(int(T1 / 60), 1) == 0):
        print(rtc.datetime)
        fz.LED(1, npx, [0, 0, 0], [0, 0, npx.v], 1)

# 2. Measurement
# 2-1. Initizalize for each measurement
    time_now, ts, Fil_Log = fz.ini_log(T1, rtc, npx)
    fz.prt(Fil_Log, "\n***********************************************************")
    fz.prt(Fil_Log, "\n\t\tFZ Measurement Log")
    fz.prt(Fil_Log, "\n***********************************************************")
    fz.prt(Fil_Log, "\n Measurement started on {} (yyyymmdd_HHMM)".format(time_now))
    fz.prt(Fil_Log, "\n\tLog file is {}".format(Fil_Log))
    fz.prt(Fil_Log, "\n\tIMU data file is IMU_BNO055_{0}.txt".format(time_now))

    if gps_i == 1 :
        fz.prt(Fil_Log, "\nGPS data ...\t{:6.1f}[s]".format(t.monotonic()-ts))
        nmea1, nmea2, lat, lon = fz.log_gps(i2c, rtc, npx)
        fz.prt(Fil_Log, "\n\tlatitude = {:9.4f}".format(lat))
        fz.prt(Fil_Log, "\n\tlongitude = {:9.4f}".format(lon))
        fz.prt(Fil_Log, "\n\tNMEA output (RMC and GGA)")
        fz.prt(Fil_Log, "\n\t {}\n".format(nmea1))
        fz.prt(Fil_Log, "\n\t {}\n".format(nmea2))
# 2-3. Logging BMP280 temperature and barometric pressure(only if sen_i == 1)
    if sen_i == 1 :
        temp, pres, hmd, Slev = fz.log_sen(i2c)
        fz.prt(Fil_Log, "\nAir data ...\t {:6.1f}[s]".format(t.monotonic()-ts))
        fz.prt(Fil_Log, "\n\tAir Temp.  : {:.2f} C".format(temp))
        fz.prt(Fil_Log, "\n\tAir Pres.  : {:.2f} hPa".format(pres))
        fz.prt(Fil_Log, "\n\tAir Humi.  : {:.2f} %".format(hmd))
        fz.prt(Fil_Log, "\n\tSound level    : {} ".format(Slev))
# 2-4. Logging IMU (the sensor can be switched by imu_i )
    if imu_i == 1 :
        fz.prt(Fil_Log, "\n\nLogging IMU ...\t{:6.1f} [s]".format(t.monotonic()-ts))
        Fil_low, Nl, imu_cal = fz.log_imu(time_now, T2, Hz1, Hz2, i2c, npx)
        fz.prt(Fil_Log, "\n\tLogging finished ...")
        fz.prt(Fil_Log, "\n\t\t The number of lines : {:5d}".format(Nl))
        fz.prt(Fil_Log, "\n\t\t ---> {}\n".format(Fil_low))
        fz.prt(Fil_Log, "\n\t Calibration / offset values ...")
        fz.prt(Fil_Log, "\n\t\t Accerelerometer offset, {} ".format(imu_cal[0:3]))
        fz.prt(Fil_Log, "\n\t\t Gyroscope       offset, {} ".format(imu_cal[3:6]))
        fz.prt(Fil_Log, "\n\t\t Magnetometer    offset, {} ".format(imu_cal[6:9]))
        fz.prt(Fil_Log, "\n\t\t Accerelerometer radius, {} ".format(imu_cal[9]))
        fz.prt(Fil_Log, "\n\t\t Magnetometer    radius, {} ".format(imu_cal[10]))
# 3. Calculate wave statistics from IMU data
    if psd_i == 1 :
        fz.prt(Fil_Log, "\nSpec. Analysis ...\t{:6.1f}[s]".format(t.monotonic()-ts))
        fz.prt(Fil_Log, "\n\tEstimates PSD ...")
        # gc.collect()
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
        rb_sent, Nr = fz.iri(Hs, Fp, Pzz, lat, lon, npx, Fil_Log)
        fz.prt(Fil_Log, "\n\tTrasfrer data =\t{}".format(rb_sent))
        fz.prt(Fil_Log, "\n\tNumber of Retry = \t{}".format(Nr))

# 5. Finished
    fz.prt(Fil_Log, "\n\nFinished after {:6.1f} [s]".format(t.monotonic()-ts))
    fz.prt(Fil_Log, "\n**********************************************************")

# 6. Reload
    storage.umount("/sd")
except OSError :
    pass
wdt.feed()
time_alarm = alarm.time.TimeAlarm(monotonic_time=t.monotonic() + 10)
alarm.exit_and_deep_sleep_until_alarms(time_alarm)
