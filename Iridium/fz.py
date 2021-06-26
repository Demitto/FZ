"""
Last updated : Jun. 7PROGRAM OF FZ, Written by T.K.       Last updated : Jun. 7, 2021
"""
import time
import board
import busio
import os
import neopixel
import ulab
from ulab import numerical as num
from ulab.vector import sin
from ulab.vector import cos
import storage
import digitalio
import adafruit_sdcard
import adafruit_pcf8523

def prt(Fil_Log, Str):
    # Objective : Make logs of the microcontroller behavior
    with open(Fil_Log, "a") as flog:
        flog.write(Str)
        flog.flush()
    # Uncomment the following line for debugging
    # print(Str)
    return

def npx_blk2(dt, pixels):
    # Objective : Blink neopixel every dt seconds with the color specified by "pixels"
    for j in range(1) :
        pixels[0] = (0, 0, 255)
        time.sleep(dt / 2)
        pixels[0] = (0, 0, 0)
        time.sleep(dt / 2)
    return pixels

def npx_blk(dt, pixels, pix_val):
    # Objective : Blink neopixel every dt seconds with the color specified by "pixels"
    for j in range(10) :
        pixels[0] = (0, pix_val, 0)
        time.sleep(dt / 2)
        pixels[0] = (0, 0, pix_val)
        time.sleep(dt / 2)
    return pixels

def rtc_now(i2c):
    # Create the RTC instance:
    rtc = adafruit_pcf8523.PCF8523(i2c)
    return rtc.datetime

def ior_prm():
    # Objective : read parameters N,dt, and T from "/param.txt"
    # N  : the number of calculation samples
    # dt : time step for data logging [s]
    # T  : total time for data logging [s]
    with open("param.txt", "r") as fp:
        for line in fp:
            if "dt" in line:
                dt = float(line[:-1].split("=")[1])
            if "T" in line:
                T = float(line[:-1].split("=")[1])
            if "Nw" in line:
                Nw = int(line[:-1].split("=")[1])
            if "Ncut" in line:
                Ncut = int(line[:-1].split("=")[1])
            if "imu" in line:
                imu_i = int(line[:-1].split("=")[1])
            if "gps" in line:
                gps_i = int(line[:-1].split("=")[1])
            if "psd" in line:
                psd_i = int(line[:-1].split("=")[1])
            if "pix_val" in line:
                pix_val = int(line[:-1].split("=")[1])
    frq = ulab.linspace(2 * 3.14 / dt / Nw, 2 * 3.14 / dt, Nw)
    return T, dt, Nw, Ncut, frq, imu_i, gps_i, psd_i, pix_val

def dvc_ini(imu_i, gps_i, set_rtc, set_time):
    # Objective : Initialization : NeoPixel, I2C, SDcard, IMU, GPS, Output Folder
    # Setup NeoPixel
    pixels = neopixel.NeoPixel(board.NEOPIXEL, 1)
    # Blink fast to show the system is initializing
    pixels = npx_blk(0.05, pixels, 255)
    # Setup I2C
    i2c = busio.I2C(board.SCL, board.SDA)
    # Setup RTC
    if(set_rtc) :
        rtc = adafruit_pcf8523.PCF8523(i2c)
        rtc.datetime = time.struct_time(set_time)
        t = rtc.datetime
        print("The date is {}/{}/{}".format(t.tm_mday, t.tm_mon, t.tm_year))
        print("The time is {}:{:02}:{:02}".format(t.tm_hour, t.tm_min, t.tm_sec))
    # Setup SD card / sdioio or adafruit_sdcard
    SD_CS = board.D10
    spi = busio.SPI(board.SCK, board.MOSI, board.MISO)
    cs = digitalio.DigitalInOut(SD_CS)
    sdcard = adafruit_sdcard.SDCard(spi, cs)
    # RTC Time Now
    # rtc_t = rtc_now(i2c)
    # time_now = "{:04d}{:02d}{:02d}_{:02d}{:02d}".format(
    #    rtc_t.tm_year, rtc_t.tm_mon, rtc_t.tm_mday, rtc_t.tm_hour, rtc_t.tm_min)
    Dir_Out = "data"
    vfs = storage.VfsFat(sdcard)
    storage.mount(vfs, "/sd")

    # Setup IMU
    if imu_i == 1 :
        import adafruit_bno055
        imu = adafruit_bno055.BNO055_I2C(i2c)
        while not imu.calibrated:
            if(num.sum(imu.calibration_status[1:4]) >= 9):
                pixels[0] = (255, 255, 255)  # NeoPixiel Orange
            elif(num.sum(imu.calibration_status[1:4]) >= 6):
                pixels[0] = (255, 255, 0)  # NeoPixiel Orange
            elif(num.sum(imu.calibration_status[1:4]) >= 3):
                pixels[0] = (125, 0, 0)  # NeoPixiel Orange
            else :
                pixels[0] = (255, 0, 0)  # NeoPixiel Orange
            time.sleep(1)
    elif imu_i == 2 :
        from adafruit_bno08x.i2c import BNO08X_I2C
        from adafruit_bno08x import (
            BNO_REPORT_ACCELEROMETER,
            BNO_REPORT_GYROSCOPE,
            BNO_REPORT_MAGNETOMETER,
            BNO_REPORT_ROTATION_VECTOR,
            BNO_REPORT_LINEAR_ACCELERATION)
        imu = BNO08X_I2C(i2c)
        imu.enable_feature(BNO_REPORT_ACCELEROMETER)
        imu.enable_feature(BNO_REPORT_GYROSCOPE)
        imu.enable_feature(BNO_REPORT_MAGNETOMETER)
        imu.enable_feature(BNO_REPORT_ROTATION_VECTOR)
        imu.enable_feature(BNO_REPORT_LINEAR_ACCELERATION)
        imu.begin_calibration()
        while not (imu.calibration_status == 3):
            if(imu.calibration_status == 2):
                pixels[0] = (255, 255, 255)  # NeoPixiel Orange
            elif(imu.calibration_status == 1):
                pixels[0] = (255, 255, 0)  # NeoPixiel Orange
            elif(imu.calibration_status == 0):
                pixels[0] = (255, 0, 0)  # NeoPixiel Orange
            time.sleep(1)
        imu.save_calibration_data()
    elif imu_i == 3 :
        import adafruit_lsm6ds.lsm6ds33
        imu = adafruit_lsm6ds.lsm6ds33.LSM6DS33(i2c)
    elif imu_i == 4 :
        from adafruit_lsm6ds.ism330dhcx import ISM330DHCX
        imu = ISM330DHCX(i2c)
    # Create Output folder
    try:
        os.mkdir("/sd/" + Dir_Out)
    except OSError:
        pass

    # Set up GPS if used
    if gps_i == 1 :
        gps_log = True
    else :
        gps_log = False

    if gps_log:
        import adafruit_gps
        # uart = busio.UART(board.TX, board.RX, baudrate=9600, timeout=10)
        # gps = adafruit_gps.GPS(uart, debug=False)
        gps = adafruit_gps.GPS_GtopI2C(i2c, debug=False)
        gps.send_command(b"PMTK314,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0")
        gps.send_command(b"PMTK220,1000")
        gps.update()
        pixels[0] = (255, 0, 0)  # NeoPixiel Red1
        while not gps.has_fix:
            pixels[0] = (255, 0, 0)  # NeoPixiel Red1
            gps.update()
            pixels[0] = (0, 0, 0)  # NeoPixiel Red1
        Dir_Out = Dir_Out + "_GPS"
    else:
        gps = 1
    # Create Output folder
    try:
        os.mkdir("/sd/" + Dir_Out)
    except OSError:
        pass
    return pixels, i2c, imu, gps, Dir_Out

def log_ini(i2c, pixels, pix_val):
    npx_blk(.1, pixels, pix_val)
    rtc_t = rtc_now(i2c)
    time_now = "{:04d}{:02d}{:02d}_{:02d}{:02d}".format(
        rtc_t.tm_year, rtc_t.tm_mon, rtc_t.tm_mday, rtc_t.tm_hour, rtc_t.tm_min)
    ts = time.monotonic()
    Fil_Log = "/sd/data/Log_" + time_now + ".txt"
    return time_now, ts, Fil_Log

def iow_imu(Dir_Out, time_now, i2c, imu, imu_i, pixels):
    # Objective : Logging IMU data every dt seconds to Fil_raw
    # Read parameters
    T, dt, Nw, Ncut, frq, imu_i, gps_i, psd_i, pix_val = ior_prm()
    # Logging duration is just 1 min less than T
    T_log = T - 120
    # Logging ...
    if (imu_i == 1) :
        Fil_raw = "/sd/{0}/IMU_BNO055_{1}.txt".format(Dir_Out, time_now)
        with open(Fil_raw, "w") as out_f:
            t_s = time.monotonic()  # Time after the feather board is on
            t_n = time.monotonic() - t_s  # Time after the logging is started
            i_l = 1  # index of the logging for checking error
            out_f.write(
                "time\t"
                "bn5_a_x\tbn5_a_y\tbn5_a_z\t"
                "bn5_gyx\tbn5_gyy\tbn5_gyz\t"
                "bn5_e_h\tbn5_e_r\tbn5_e_p\t"
                "bn5_g_x\tbn5_g_y\tbn5_g_z\t"
                "cal_1\tcal_2\tcal_3\tcal_4\r\n"
            )
            out_f.flush()
            while t_n < T_log + dt:
                t_n = time.monotonic() - t_s
                if t_n >= dt * float(i_l):
                    t_n = time.monotonic() - t_s
                    out_f.write("%10.5f\t" % t_n)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.linear_acceleration)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.gyro)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.euler)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.gravity)
                    out_f.write("%1d\t%1d\t%1d\t%1d\r\n" % imu.calibration_status)
                    i_l += 1
            out_f.flush()
    elif (imu_i == 2) :
        Fil_raw = "/sd/{0}/IMU_BNO08x_{1}.txt".format(Dir_Out, time_now)
        # Logging ...
        with open(Fil_raw, "w") as out_f:
            t_s = time.monotonic()  # Time after the feather board is on
            t_n = time.monotonic() - t_s  # Time after the logging is started
            i_l = 1  # index of the logging for checking error
            out_f.write(
                "time\t"
                "bn8_a_x\tbn8_a_y\tbn8_a_z\t"
                "bn8_l_x\tbn8_l_y\tbn8_l_z\t"
                "bn8_gyx\tbn8_gyy\tbn8_gyz\tbn8cali\r\n"
            )
            out_f.flush()
            while t_n < T_log + dt:
                t_n = time.monotonic() - t_s
                if t_n >= dt * float(i_l):
                    t_n = time.monotonic() - t_s
                    out_f.write("%10.5f\t" % t_n)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.acceleration)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.linear_acceleration)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.gyro)
                    out_f.write("%1d\r\n" % imu.calibration_status)
                    i_l += 1
            out_f.flush()
    elif (imu_i == 3) :
        Fil_raw = "/sd/{0}/IMU_LSM6DS33_{1}.txt".format(Dir_Out, time_now)
        with open(Fil_raw, "w") as out_f:
            t_s = time.monotonic()  # Time after the feather board is on
            t_n = time.monotonic() - t_s  # Time after the logging is started
            i_l = 1  # index of the logging for checking error
            out_f.write(
                "time\t"
                "lsm_a_x\tlsm_a_y\tlsm_a_z\t"
                "lsm_gyx\tlsm_gyy\tlsm_gyz\r\n"
            )
            out_f.flush()
            while t_n < T_log + dt:
                t_n = time.monotonic() - t_s
                if t_n >= dt * float(i_l):
                    t_n = time.monotonic() - t_s
                    out_f.write("%10.5f\t" % t_n)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.acceleration)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\r\n" % imu.gyro)
                    i_l += 1
            out_f.flush()
    elif (imu_i == 4) :
        Fil_raw = "/sd/{0}/IMU_ISM330DXCX_{1}.txt".format(Dir_Out, time_now)
        with open(Fil_raw, "w") as out_f:
            t_s = time.monotonic()  # Time after the feather board is on
            t_n = time.monotonic() - t_s  # Time after the logging is started
            i_l = 1  # index of the logging for checking error
            out_f.write(
                "time\t"
                "ism_a_x\tism_a_y\tism_a_z\t"
                "ism_gyx\tism_gyy\tism_gyz\r\n"
            )
            out_f.flush()
            while t_n < T_log + dt:
                t_n = time.monotonic() - t_s
                if t_n >= dt * float(i_l):
                    t_n = time.monotonic() - t_s
                    out_f.write("%10.4f\t" % t_n)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\t" % imu.acceleration)
                    out_f.write("%8.4f\t%8.4f\t%8.4f\r\n" % imu.gyro)
                    i_l += 1
            out_f.flush()
    N_l = lin_nmb(Fil_raw)
    return Fil_raw, time_now, N_l

def flt_low(Fil_raw, N_l, dt):
    # only for bno055 yet
    Fil_low = Fil_raw[0:-4]+'_low.txt'
    N2 = int(1/dt/2)
    N1 = int(N_l/N2)
    t1 = ulab.zeros(N2)
    ax, ay, az = ulab.zeros(N2), ulab.zeros(N2), ulab.zeros(N2)
    gx, gy, gz = ulab.zeros(N2), ulab.zeros(N2), ulab.zeros(N2)
    hd, rl, pc = ulab.zeros(N2), ulab.zeros(N2), ulab.zeros(N2)
    with open(Fil_raw, "r") as in_f:
        header = in_f.readline()
        with open(Fil_low, "w") as out_f:
            out_f.write(header)
            for i in range(N1) :
                for j in range(N2) :
                    dat = in_f.readline().split("\t")
                    t1[j] = float(dat[0])
                    ax[j], ay[j], az[j] = float(dat[1]), float(dat[2]), float(dat[3])
                    gx[j], gy[j], gz[j] = float(dat[4]), float(dat[5]), float(dat[6])
                    hd[j], rl[j], pc[j] = float(dat[7]), float(dat[8]), float(dat[9])
                out_f.write("%10.5f\t" % num.mean(t1))
                out_f.write("%8.4f\t" % num.mean(ax))
                out_f.write("%8.4f\t" % num.mean(ay))
                out_f.write("%8.4f\t" % num.mean(az))
                out_f.write("%8.4f\t" % num.mean(gx))
                out_f.write("%8.4f\t" % num.mean(gy))
                out_f.write("%8.4f\t" % num.mean(gz))
                out_f.write("%8.4f\t" % num.median(hd))
                out_f.write("%8.4f\t" % num.mean(rl))
                out_f.write("%8.4f\r\n" % num.mean(pc))
    return Fil_low, N1

def lin_nmb(F_IMU):
    # Objective : Calculate the number of lines in the file "F_IMU"
    file = open(F_IMU, "r")
    N_l = 0
    for line in file:
        if line != "\n":
            N_l += 1
    file.close()
    return N_l

def psd_seg(F_IMU, N_low):
    # Objective : Read raw data and calculate the psd for "Nseg" segment"s".
    # The window size is N. The coordinate is converted to Earth Coordinate
    # following Bender et al., 2010.
    T, dt, Nw, Ncut, frq, imu_i, gps_i, psd_i, pix_val = ior_prm()
    d2r = 3.1415 / 180.0
    ax, ay, az = ulab.zeros(Nw), ulab.zeros(Nw), ulab.zeros(Nw)
    hd, rl, pc = ulab.zeros(Nw), ulab.zeros(Nw), ulab.zeros(Nw)
    d_ax, d_ay, d_az = ulab.zeros(Nw), ulab.zeros(Nw), ulab.zeros(Nw)
    Pxx, Pyy, Pzz = ulab.zeros(Nw), ulab.zeros(Nw), ulab.zeros(Nw)
    Qxz, Qyz, Cxy = ulab.zeros(Nw), ulab.zeros(Nw), ulab.zeros(Nw)
    with open(F_IMU, "r") as file:
        Nseg = ulab.vector.floor((N_low - 1) / Nw)
        next(file)
        for seg in range(Nseg):
            for i in range(Nw):
                dat = file.readline().split("\t")
                d_ax[i], d_ay[i], d_az[i] = float(dat[1]), float(dat[2]), float(dat[3])
                hd[i], rl[i], pc[i] = float(dat[7]), float(dat[8]), float(dat[9])
            hd, rl, pc = hd * d2r, rl * d2r, pc * d2r
            ax = (d_ax * cos(pc) * cos(hd) +
                  d_ay * (sin(rl) * sin(pc) * cos(hd) - cos(rl) * sin(hd)) +
                  d_az * (cos(rl) * sin(pc) * cos(hd) + sin(rl) * sin(hd)))
            ay = (d_ax * cos(pc) * sin(hd) +
                  d_ay * (sin(rl) * sin(pc) * sin(hd) + cos(rl) * cos(hd)) +
                  d_az * (cos(rl) * sin(pc) * sin(hd) - sin(rl) * cos(hd)))
            az = -d_ax * sin(pc) + d_ay * sin(rl) * cos(pc) + d_az * cos(rl) * cos(pc)
            Pxx += psd(frq, ax) / float(Nseg)
            Pyy += psd(frq, ay) / float(Nseg)
            Pzz += psd(frq, az) / float(Nseg)
            Qxz += qsd(frq, ax, az, 2) / float(Nseg)
            Qyz += qsd(frq, ay, az, 2) / float(Nseg)
            Cxy += qsd(frq, ax, ay, 1) / float(Nseg)

    Fil_psd = F_IMU[0:-4] + "_psd.txt"
    with open(Fil_psd, "w") as f_psd:
        for j in range(Nw/2) :
            f_psd.write("%7.3f\t %7.3f\r\n" % (frq[j], Pzz[j]))
        f_psd.flush()

    return Pxx, Pyy, Pzz, Qxz, Qyz, Cxy

def psd(frq, X):
    # Objective : Calculate the power spectral density (psd)
    ffr, ffi = ulab.zeros(len(X)), ulab.zeros(len(X))
    X = X - num.mean(X)
    ffr, ffi = ulab.fft.fft(X)
    ffr, ffi = ffr / len(X), ffi / len(X)
    ffr, ffi = ffr / (frq ** 2), ffi / (frq ** 2)
    # Save the memory to use "ffr" for the output
    ffr = ffr ** 2 + ffi ** 2
    ffr[int(len(X) / 2) :] = 0
    ffr = 2 * ffr
    return ffr

def qsd(frq, X, Y, flag):
    # Objective : Calculate the psd and quadrature spectral density (qsd).
    ffrx, ffix = ulab.zeros(len(X)), ulab.zeros(len(X))
    ffry, ffiy = ulab.zeros(len(X)), ulab.zeros(len(X))
    A = ulab.zeros(len(X))
    X = X - num.mean(X)
    Y = Y - num.mean(Y)
    ffrx, ffix = ulab.fft.fft(X)
    ffrx, ffix = ffrx / len(X), ffix / len(X)
    ffrx, ffix = ffrx / (frq ** 2), ffix / (frq ** 2)
    ffry, ffiy = ulab.fft.fft(Y)
    ffry, ffiy = ffry / len(Y), ffiy / len(Y)
    ffry, ffiy = ffry / (frq ** 2), ffiy / (frq ** 2)
    # Calculate the Quadrature
    # qst+i*qsi = ( xr - i*xi ) * (yr + i*yr)
    #           =  ( xr*yr + xi*yi) + i (-xi*yr + xr*yi)
    if flag == 1:
        A = 2 * ( ffrx * ffry + ffix * ffiy )
    elif flag == 2:
        A = 2 * ( ffrx * ffiy - ffix * ffry )
    A[int(len(A)/2):] = 0
    return A

def blk_sta(Pxx, Pyy, Pzz, Qxz, Qyz, Cxy):
    # Objective : Calculate the bulk wave statistics based on the psd & csd
    T, dt, Nw, Ncut, frq, imu_i, gps_i, psd_i, pix_val = ior_prm()
    a1, b1 = ulab.zeros(len(Pxx)), ulab.zeros(len(Pxx))
    Mwd, Mds = ulab.zeros(len(Pxx)), ulab.zeros(len(Pxx))
    Hs = 4 * (num.sum(Pzz[Ncut:]) ** 0.5)
    Hsx = 4 * (num.sum(Pxx[Ncut:]) ** 0.5)
    Hsy = 4 * (num.sum(Pyy[Ncut:]) ** 0.5)
    Fp = frq[Ncut-1+ulab.numerical.argmax(Pzz[Ncut:])] / 6.28
    Fpx = frq[Ncut-1+ulab.numerical.argmax(Pxx[Ncut:])] / 6.28
    Fpy = frq[Ncut-1+ulab.numerical.argmax(Pyy[Ncut:])] / 6.28
    a1 = Qxz / (((Pxx + Pyy) * Pzz) ** 0.5)
    b1 = Qyz / (((Pxx + Pyy) * Pzz) ** 0.5)
    # a2 = (Pxx - Pyy) / (Pxx + Pyy)
    # b2 = (2 * Cxy) / (Pxx + Pyy)
    Mwd = ulab.vector.atan(b1 / a1)
    Mds = (2 * (1 - (a1 ** 2 + b1 ** 2) ** 0.5)) ** 0.5
    return Hs, Hsx, Hsy, Fp, Fpx, Fpy, Mwd, Mds

def iri(data, length)
    
