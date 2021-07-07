"""
Last updated : Jun. 7PROGRAM OF FZ, Written by T.K.       Last updated : Jul. 5, 2021
"""
import time
import board
import busio
import os
import neopixel
import array
import ulab
from ulab import numerical as num
from ulab.vector import sin
from ulab.vector import cos
import storage
import digitalio
import adafruit_sdcard
import adafruit_pcf8523
from micropython import const

def prt(Fil_Log, Str):
    # Objective : Make logs of the microcontroller behavior
    with open(Fil_Log, "a") as flog:
        flog.write(Str)
        flog.flush()
    # Uncomment the following line for debugging
    # print(Str)
    return

def npx_blk(dt, pix, col1, col2):
    # Objective : Blink neopixel every dt seconds with the color specified by "pix"
    for j in range(5) :
        pix[0] = col1
        time.sleep(dt / 2)
        pix[0] = col2
        time.sleep(dt / 2)
    return pix

def ior_prm():
    # Objective : read input parameters from "/param.txt"
    with open("param.txt", "r") as fp:
        for line in fp:
            if "T_imu" in line:
                T_imu = float(line[:-1].split("=")[1])
            elif "T" in line:
                T = float(line[:-1].split("=")[1])
            if "Hz1" in line:
                Hz1 = int(line[:-1].split("=")[1])
            if "Hz2" in line:
                Hz2 = int(line[:-1].split("=")[1])
            if "Nw" in line:
                Nw = int(line[:-1].split("=")[1])
            if "Ncut" in line:
                Ncut = int(line[:-1].split("=")[1])
            if "imu" in line:
                imu_i = int(line[:-1].split("=")[1])
            if "cal" in line:
                cal_i = int(line[:-1].split("=")[1])
            if "gps" in line:
                gps_i = int(line[:-1].split("=")[1])
            if "sen" in line:
                sen_i = int(line[:-1].split("=")[1])
            if "psd" in line:
                psd_i = int(line[:-1].split("=")[1])
            if "pix_val" in line:
                pval = int(line[:-1].split("=")[1])
    frq = ulab.linspace(2 * 3.14 * Hz2 / Nw, 2 * 3.14 * Hz2, Nw)
    return T, T_imu, Hz1, Hz2, Nw, Ncut, frq, imu_i, cal_i, gps_i, sen_i, psd_i, pval

def dvc_ini(T, imu_i, cal_i, gps_i, sen_i, pval, set_rtc, set_time):
    # Objective : Initialization : NeoPixel, I2C, SDcard, Output directory, Sensors

    # Setup NeoPixel--->Blink fast to show the system is initializing
    pix = neopixel.NeoPixel(board.NEOPIXEL, 1)
    pix = npx_blk(0.1, pix, [0, pval, 0], [0, 0, pval])

    # Setup I2C
    i2c = busio.I2C(board.SCL, board.SDA)

    # Setup RTC
    rtc = adafruit_pcf8523.PCF8523(i2c)
    if(set_rtc) :
        rtc.datetime = time.struct_time(set_time)

    # Setup SD card / sdioio or adafruit_sdcard
    SD_CS = board.D10
    spi = busio.SPI(board.SCK, board.MOSI, board.MISO)
    cs = digitalio.DigitalInOut(SD_CS)
    sdcard = adafruit_sdcard.SDCard(spi, cs)
    Dir_Out = "data"
    vfs = storage.VfsFat(sdcard)
    storage.mount(vfs, "/sd")

    # Create Output folder
    try:
        os.mkdir("/sd/" + Dir_Out)
    except OSError:
        pass

    # Setup IMU
    if imu_i == 1 :
        import adafruit_bno055
        imu = adafruit_bno055.BNO055_I2C(i2c)
        imu.mode = adafruit_bno055.NDOF_FMC_OFF_MODE
        if cal_i == 1 :
            while not imu.calibrated:
                if(num.sum(imu.calibration_status[1:4]) >= 9):
                    pix[0] = (pval, pval, pval)  # NeoPixiel Orange
                elif(num.sum(imu.calibration_status[1:4]) >= 6):
                    pix[0] = (pval, pval, 0)  # NeoPixiel Orange
                elif(num.sum(imu.calibration_status[1:4]) >= 3):
                    pix[0] = (pval/2, 0, 0)  # NeoPixiel Orange
                else :
                    pix[0] = (pval, 0, 0)  # NeoPixiel Orange
                time.sleep(.2)
                print(imu.calibration_status)
            pix[0] = (0, 0, pval)  # NeoPixiel Orange
    else:
        imu = []

    # Setup other sensors if used
    if sen_i == 1 :
        import adafruit_bmp280
        bmp = adafruit_bmp280.Adafruit_BMP280_I2C(i2c)
        import audiobusio
        mic = audiobusio.PDMIn(board.MICROPHONE_CLOCK, board.MICROPHONE_DATA,
                               sample_rate=16000, bit_depth=16)
        micv = array.array('H', [0] * 160)
        mic.record(micv, len(micv))
        import adafruit_sht31d
        hmd = adafruit_sht31d.SHT31D(i2c)
    else:
        bmp, mic, micv, hmd = [], [], [], []

    # Set up GPS if used
    pix[0] = [pval, 0, 0]  # NeoPixiel Red1
    if gps_i == 1 :
        import adafruit_gps
        gps = adafruit_gps.GPS_GtopI2C(i2c, debug=False)
        pass_i = 0
        while (pass_i == 0) :
            try :
                gps.send_command(b'PMTK225,0')
                pass_i = 1
            except OSError:
                time.sleep(.45)
        pix[0] = (pval, int(pval/4), 0)  # NeoPixiel Red1
        while not (gps.has_fix) :
            gps.update()
        gps.send_command(b'PMTK314,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        gps.send_command(b"PMTK220,1000")
        for i in range(4) :
            print(gps.readline())
        rtc.datetime = time.struct_time(gps.timestamp_utc)
    else:
        gps = []
    return pix, i2c, rtc, imu, gps, bmp, mic, micv, hmd

def gps_log(i2c, gps, rtc, pix, pval) :
    pix[0] = (pval, int(pval/4), 0)  # NeoPixiel Red1
    pass_i = 0
    while (pass_i == 0) :
        try :
            gps.send_command(b"PMTK225,0")
            gps.update()
            pass_i = 1
        except OSError:
            time.sleep(.45)
    while not (gps.has_fix) :
        gps.update()
    for i in range(10) :
        gps.update()
        nmea = str(gps.readline(), "ascii").strip()
        # print(gps.update(), gps.has_fix, nmea, gps.latitude, gps.longitude)
    # print(gps.datetime, '\r\n')
    lat = gps.latitude
    lon = gps.longitude
    gps.send_command(b"PMTK225,1,60000,1140000,360000000,300000")
    pix[0] = (0, 0, pval)  # NeoPixiel Blue
    return nmea, lat, lon

def log_ini(i2c, rtc, pix, pval):
    npx_blk(.5, pix, [0, pval, 0], [0, 0, pval])
    time_now = "{:04d}{:02d}{:02d}_{:02d}{:02d}{:02d}".format(
        rtc.datetime.tm_year, rtc.datetime.tm_mon, rtc.datetime.tm_mday,
        rtc.datetime.tm_hour, rtc.datetime.tm_min, rtc.datetime.tm_sec)
    ts = time.monotonic()
    Fil_Log = "/sd/data/Log_" + time_now + ".txt"
    return time_now, ts, Fil_Log

def iow_imu(time_now, T, T_imu, Hz1, Hz2, i2c, imu, imu_i, pix):
    # Objective : Logging IMU data every dt seconds to Fil_low
    Dir_Out = "data"
    if (imu_i == 1) :
        Fil_low = "/sd/{0}/IMU_BNO055_{1}.txt".format(Dir_Out, time_now)
        with open(Fil_low, "w") as out_f:
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
            i, tr, acc, gyr, eul, grv, cal = 0, 0, [], [], [], [], []
            while t_n < T_imu + 1/Hz1 :
                t_n = time.monotonic() - t_s
                if t_n >= float(i_l)/Hz1 :
                    i_l += 1
                    i += 1
                    tr += t_n
                    acc += imu.linear_acceleration
                    gyr += imu.gyro
                    eul += imu.euler
                    grv += imu.gravity
                    cal += imu.calibration_status
                    if((i_l / Hz1 * Hz2) % 1 == 0) :
                        tr = tr/i
                        # print([i, i_l, tr])
                        acc = [sum(acc[j::3])/i for j in range(3)]
                        gyr = [sum(gyr[j::3])/i for j in range(3)]
                        eul = [sum(eul[j::3])/i for j in range(3)]
                        grv = [sum(grv[j::3])/i for j in range(3)]
                        cal = [sum(cal[j::3])/i for j in range(4)]
                        # print(acc)
                        out_f.write("%10.5f\t" % tr)
                        out_f.write("%8.4f\t%8.4f\t%8.4f\t" % tuple(acc))
                        out_f.write("%8.4f\t%8.4f\t%8.4f\t" % tuple(gyr))
                        out_f.write("%8.4f\t%8.4f\t%8.4f\t" % tuple(eul))
                        out_f.write("%8.4f\t%8.4f\t%8.4f\t" % tuple(grv))
                        out_f.write("%1d\t%1d\t%1d\t%1d\r\n" % tuple(cal))
                        i, tr, acc, gyr, eul, grv, cal = 0, 0, [], [], [], [], []
            out_f.flush()
        file = open(Fil_low, "r")
        N_l = 0
        for line in file:
            if line != "\n":
                N_l += 1
        file.close()
    return Fil_low, time_now, N_l

def imu_chg(imu, mode_i):
    if mode_i == 0 :
        # imu.accel_mode = adafruit_bno055.ACCEL_SUSPEND_MODE
        # imu.gyro_mode = adafruit_bno055.GYRO_SUSPEND_MODE
        # imu.magnet_mode = adafruit_bno055.MAGNET_SUSPEND_MODE
        imu._write_register(const(0x3E), const(0x02))
    elif mode_i == 1 :
        # imu.accel_mode = adafruit_bno055.ACCEL_NORMAL_MODE
        # imu.gyro_mode = adafruit_bno055.GYRO_NORMAL_MODE
        # imu.magnet_mode = adafruit_bno055.MAGNET_REGULAR_MODE
        imu._write_register(const(0x3E), const(0x00))


def psd_seg(F_IMU, N_low, Nw, frq) :
    # Objective : Read raw data and calculate the psd for "Nseg" segment"s".
    # The window size is N. The coordinate is converted to Earth Coordinate
    # following Bender et al., 2010.
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
        A = 2 * (ffrx * ffry + ffix * ffiy)
    elif flag == 2:
        A = 2 * (ffrx * ffiy - ffix * ffry)
    A[int(len(A)/2):] = 0
    return A

def blk_sta(Ncut, frq, Pxx, Pyy, Pzz, Qxz, Qyz, Cxy):
    # Objective : Calculate the bulk wave statistics based on the psd & csd
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
