""" ############################################################################
#    Subroutine programs OF FZ, Written by T.K.      Last updated : Nov. 2, 2021
#
# -Overbies of each subroutines
#   -prt(Fil_Log, Str) :
#       - Logging the Str (MCU actions) to Fil_Log (serial cosoal if uncommented).
#   -LED(dt, npx, col1, col2, Nb) :
#       - Blink the neopixel. The interval dt and color col1/col2 and the iteration Nb
#   -prm() :
#       - Install the parameter settings from param.txt
#   -ini_all(imu_i, cal_i, gps_i, sen_i, set_rtc, set_time) :
#       - Initialize all the selected devices based on ???_i index from prm.
#   -ini_log(T1, imu, rtc, npx) :
#       - Initializaion for each measurements (waiting, gps, etc)
#   -log_gps(gps, rtc, npx) :
#       - Logging gps (only at once in each measurement cycles).
#   -log_imu(time_now, T2, Hz1, Hz2, imu, mgn, npx) :
#       - Logging imu at Hz2 (sampled at Hz1 and subsample Hz2 via ave.)
#   -cal_spc(F_IMU, Nl, Nw, frq) :
#      -Spectral Analysis based on the following three types of the calculation.
#      -cal_psd(frq, X):
#           - Calculate power spectral density
#      -cal_qsd(frq, X, Y, flag):
#           - Calculate quadrature spectral density
#      -cal_blk(Ncut, frq, Pxx, Pyy, Pzz, Qxz, Qyz, Cxy):
#           - Calculate bulk wave statistics
#   -iri(rb, Hs, Fp, Pzz, lat, lon, npx, Fil_Log) :
#       - Satellite Data Transfer using Iridium-SBD (RockBlock 9603)
#   ############################################################################
"""
import time
import board
import busio
import os
import neopixel
import array
from ulab import numerical as num
from ulab import zeros, linspace, vector, fft
from ulab.vector import sin, cos
import storage
import digitalio
import adafruit_sdcard
import adafruit_pcf8523
from micropython import const
import struct

def prt(Fil_Log, Str):
    with open(Fil_Log, "a") as flog:
        flog.write(Str)
        flog.flush()
    # Uncomment the following line for debugging
    # print(Str)
    return


def LED(dt, npx, col1, col2, Nb):
    for j in range(Nb):
        npx[0] = col1
        time.sleep(dt / 2)
        npx[0] = col2
        time.sleep(dt / 2)
    return npx


def prm():
    with open("param.txt", "r") as fp:
        for line in fp:
            if "T1" in line:
                T1 = float(line[:-1].split("=")[1])
            elif "T2" in line:
                T2 = float(line[:-1].split("=")[1])
            if "Hz1" in line:
                Hz1 = int(line[:-1].split("=")[1])
            if "Hz2" in line:
                Hz2 = int(line[:-1].split("=")[1])
            if "Nw" in line:
                Nw = int(line[:-1].split("=")[1])
            if "Nc" in line:
                Nc = int(line[:-1].split("=")[1])
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
            if "iri" in line:
                iri_i = int(line[:-1].split("=")[1])
    frq = linspace(1e-10, 2 * 3.14 * Hz2 * (Nw - 1) / Nw, Nw)
    return T1, T2, Hz1, Hz2, Nw, Nc, frq, imu_i, cal_i, gps_i, sen_i, psd_i, iri_i


def ini_all(imu_i, cal_i, gps_i, sen_i, iri_i, set_rtc, set_time):
    # Setup NeoPixel--->Blink fast to show the system is initializing
    npx = neopixel.NeoPixel(board.NEOPIXEL, 1)
    npx.v = 250
    npx = LED(0.2, npx, [0, 0, npx.v], [0, npx.v, 0], 10)
    # Setup I2C
    i2c = busio.I2C(board.SCL, board.SDA)
    # Setup RTC
    rtc = adafruit_pcf8523.PCF8523(i2c)
    if set_rtc:
        rtc.datetime = time.struct_time(set_time)
    # Setup SD card / sdioio or adafruit_sdcard
    SD_CS = board.D10
    spi = busio.SPI(board.SCK, board.MOSI, board.MISO)
    cs = digitalio.DigitalInOut(SD_CS)
    sdcard = adafruit_sdcard.SDCard(spi, cs)
    vfs = storage.VfsFat(sdcard)
    storage.mount(vfs, "/sd")
    # Setup IMU
    if imu_i == 1:
        import adafruit_bno055

        imu = adafruit_bno055.BNO055_I2C(i2c)
        # imu.mode = adafruit_bno055.NDOF_FMC_OFF_MODE
        if cal_i == 1:
            while not imu.calibrated:
                npx = LED(0.1, npx, [0, 0, npx.v], [npx.v / 3, npx.v / 3, npx.v / 3], 2)
                print(imu.calibration_status)
        mgn = []
    else:
        imu = []
        mgn = []
    # Setup other sensors if used
    if sen_i == 1:
        import adafruit_bmp280

        bmp = adafruit_bmp280.Adafruit_BMP280_I2C(i2c)
        import audiobusio

        mic = audiobusio.PDMIn(
            board.MICROPHONE_CLOCK,
            board.MICROPHONE_DATA,
            sample_rate=16000,
            bit_depth=16,
        )
        micv = array.array("H", [0] * 160)
        mic.record(micv, len(micv))
        import adafruit_sht31d

        hmd = adafruit_sht31d.SHT31D(i2c)
    else:
        bmp, mic, micv, hmd = [], [], [], []
    # Set up GPS if used
    if gps_i == 1:
        import adafruit_gps

        gps = adafruit_gps.GPS_GtopI2C(i2c, debug=False)
        gps.gps_switch = digitalio.DigitalInOut(board.A2)
        gps.gps_switch.direction = digitalio.Direction.OUTPUT
        nmea1, nmea2, lat, lon = log_gps(i2c, gps, rtc, npx)
    else:
        gps = []

    # Setup Lockblock iridium if used
    if iri_i == 1:
        import adafruit_rockblock

        uart = board.UART()
        uart.baudrate = 19200
        rb = adafruit_rockblock.RockBlock(uart)
        rb.switch = digitalio.DigitalInOut(board.A1)
        rb.switch.direction = digitalio.Direction.OUTPUT
        rb.switch.value = False
    else:
        rb = []

    return npx, i2c, rtc, imu, mgn, gps, bmp, mic, micv, hmd, gps_i, rb


def ini_log(T1, imu, rtc, npx):
    LED(0.5, npx, [0, 0, npx.v], [0, npx.v, 0], 5)
    time_now = "{:04d}{:02d}{:02d}_{:02d}{:02d}{:02d}".format(
        rtc.datetime.tm_year,
        rtc.datetime.tm_mon,
        rtc.datetime.tm_mday,
        rtc.datetime.tm_hour,
        rtc.datetime.tm_min,
        rtc.datetime.tm_sec,
    )
    ts = time.monotonic()
    # Create Output folder
    Dir_Out = time_now[0:8]
    try:
        os.mkdir("/sd/" + Dir_Out)
    except OSError:
        pass
    Fil_Log = "/sd/" + time_now[0:8] + "/Log_" + time_now + ".txt"
    # Sleep IMU and wait by the time of its logging
    imu._write_register(const(0x3E), const(0x02))
    while not (rtc.datetime.tm_min % max(int(T1 / 60), 1) == 0):
        LED(1, npx, [0, 0, 0], [0, 0, npx.v], 1)
    imu._write_register(const(0x3E), const(0x00))
    npx[0] = (0, npx.v, 0)
    return time_now, ts, Fil_Log


def log_gps(gps, rtc, npx):
    gps.gps_switch.value = True
    pass_i = 0
    while pass_i == 0:
        LED(1.0, npx, [0, npx.v, 0], [npx.v / 3, npx.v / 3, npx.v / 3], 1)
        try:
            gps.send_command(b"PMTK225,0")
            gps.send_command(b"PMTK314,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0")
            gps.send_command(b"PMTK220,1000")
            pass_i = 1
        except OSError:
            pass
        npx[0] = (0, npx.v, 0)  # NeoPixiel Blue
    gps_type = []
    while gps_type != "$GNRMC":
        LED(0.25, npx, [0, npx.v, 0], [npx.v / 3, npx.v / 3, npx.v / 3], 1)
        nmea1 = str(gps.readline(), "ascii").strip()
        if len(nmea1.split(",")) > 6:
            if nmea1.split(",")[2] == "A":
                gps_type = nmea1.split(",")[0]
    nmea2 = str(gps.readline(), "ascii").strip()
    if nmea1.split(",")[4] == "N":
        lat = float(nmea1.split(",")[3])
    elif nmea1.split(",")[4] == "S":
        lat = -float(nmea1.split(",")[3])
    else:
        lat = -100
    if nmea1.split(",")[6] == "E":
        lon = float(nmea1.split(",")[5])
    elif nmea1.split(",")[6] == "W":
        lon = -float(nmea1.split(",")[5])
    else:
        lon = -1000
    year = int(nmea1.split(",")[9][4:6]) + 2000
    mont = int(nmea1.split(",")[9][2:4])
    day_ = int(nmea1.split(",")[9][0:2])
    hour = int(nmea1.split(",")[1][0:2])
    mins = int(nmea1.split(",")[1][2:4])
    secs = int(nmea1.split(",")[1][4:6])
    gps_time = (year, mont, day_, hour, mins, secs, 0, -1, -1)
    rtc.datetime = time.struct_time(gps_time)
    gps.send_command(b"PMTK225,4")
    gps.gps_switch.value = False
    return nmea1, nmea2, lat, lon


def log_imu(time_now, T2, Hz1, Hz2, imu, mgn, npx):
    npx[0] = [0, 0, npx.v]
    Dir_Out = time_now[0:8]
    Fil_low = "/sd/{0}/IMU_BNO055_{1}.txt".format(Dir_Out, time_now)
    Nl = 0
    with open(Fil_low, "w") as out_f:
        t_s = time.monotonic()  # Time starting the IMU logging
        t_n = time.monotonic() - t_s  # Time from the logging start
        i_l = 1  # index of the logging for checking error
        out_f.write(
            "time\t"
            "bn5_a_x\tbn5_a_y\tbn5_a_z\t"
            "bn5_e_h\tbn5_e_r\tbn5_e_p\t"
            "bn5_N_i\t"
            "bn5_ca1\tbn5_ca2\tbn5_ca3\tbn5_ca4\r\n"
        )
        out_f.flush()
        i, acc, eul, cal = 0, [], [], []
        while t_n < T2 + 1 / Hz1:
            t_n = time.monotonic() - t_s
            if t_n >= float(i_l) / Hz1:
                i_l += 1
                i += 1
                acc += imu.linear_acceleration
                if (i_l / Hz1 * Hz2) % 1 == 0:
                    acc = [sum(acc[j::3]) / i for j in range(3)]
                    eul = imu.euler  # [sum(gyr[j::3])/i for j in range(3)]
                    cal = imu.calibration_status
                    try:
                        out_f.write("%10.5f\t" % t_n)
                        out_f.write("%8.4f\t%8.4f\t%8.4f\t" % tuple(acc))
                        out_f.write("%8.4f\t%8.4f\t%8.4f\t" % tuple(eul))
                        out_f.write("%2d\t" % i)
                        out_f.write("%1d\t%1d\t%1d\t%1d\r\n" % tuple(cal))
                        out_f.flush()
                    except OSError:
                        pass
                    # print(t_n,i)
                    i, acc, eul, cal = 0, [], [], []
        Nl = int(i_l / Hz1 * Hz2)
    return Fil_low, Nl


def cal_spc(F_IMU, Nl, Nw, frq):
    d2r = 3.1415 / 180.0
    ax, ay, az = zeros(Nw), zeros(Nw), zeros(Nw)
    hd, rl, pc = zeros(Nw), zeros(Nw), zeros(Nw)
    d_ax, d_ay, d_az = zeros(Nw), zeros(Nw), zeros(Nw)
    Pxx, Pyy, Pzz = zeros(Nw), zeros(Nw), zeros(Nw)
    Qxz, Qyz, Cxy = zeros(Nw), zeros(Nw), zeros(Nw)
    with open(F_IMU, "r") as file:
        Nseg = vector.floor((Nl - 1) / Nw)
        next(file)
        for seg in range(Nseg):
            for i in range(Nw):
                dat = file.readline().split("\t")
                d_ax[i], d_ay[i], d_az[i] = float(dat[1]), float(dat[2]), float(dat[3])
                hd[i], rl[i], pc[i] = float(dat[4]), float(dat[5]), float(dat[6])
            hd, rl, pc = hd * d2r, rl * d2r, pc * d2r
            ax = (
                d_ax * cos(pc) * cos(hd)
                + d_ay * (sin(rl) * sin(pc) * cos(hd) - cos(rl) * sin(hd))
                + d_az * (cos(rl) * sin(pc) * cos(hd) + sin(rl) * sin(hd))
            )
            ay = (
                d_ax * cos(pc) * sin(hd)
                + d_ay * (sin(rl) * sin(pc) * sin(hd) + cos(rl) * cos(hd))
                + d_az * (cos(rl) * sin(pc) * sin(hd) - sin(rl) * cos(hd))
            )
            az = -d_ax * sin(pc) + d_ay * sin(rl) * cos(pc) + d_az * cos(rl) * cos(pc)
            Pxx += cal_psd(frq, ax) / float(Nseg)
            Pyy += cal_psd(frq, ay) / float(Nseg)
            Pzz += cal_psd(frq, az) / float(Nseg)
            Qxz += cal_qsd(frq, ax, az, 2) / float(Nseg)
            Qyz += cal_qsd(frq, ay, az, 2) / float(Nseg)
            Cxy += cal_qsd(frq, ax, ay, 1) / float(Nseg)
    Fil_psd = F_IMU[0:-4] + "_psd.txt"
    with open(Fil_psd, "w") as f_psd:
        f_psd.write("Freq.[Hz]\t PSD[m^2/Hz]\r\n")
        for j in range(Nw / 2):
            f_psd.write(
                "%10.8f\t %10.8f\r\n" % (frq[j] / 6.283, Pzz[j] / frq[1] * 6.283)
            )
        f_psd.flush()
    return Pxx, Pyy, Pzz, Qxz, Qyz, Cxy


def cal_psd(frq, X):
    ffr, ffi = zeros(len(X)), zeros(len(X))
    X = X - num.mean(X)
    ffr, ffi = fft.fft(X)
    ffr, ffi = ffr / len(X), ffi / len(X)
    ffr, ffi = ffr / (frq ** 2), ffi / (frq ** 2)
    # Save the memory to use "ffr" for the output
    ffr = ffr ** 2 + ffi ** 2
    ffr[int(len(X) / 2) :] = 0
    ffr = 2 * ffr
    return ffr


def cal_qsd(frq, X, Y, flag):
    # Objective : Calculate the psd and quadrature spectral density (qsd).
    ffrx, ffix = zeros(len(X)), zeros(len(X))
    ffry, ffiy = zeros(len(X)), zeros(len(X))
    A = zeros(len(X))
    X = X - num.mean(X)
    Y = Y - num.mean(Y)
    ffrx, ffix = fft.fft(X)
    ffrx, ffix = ffrx / len(X), ffix / len(X)
    ffrx, ffix = ffrx / (frq ** 2), ffix / (frq ** 2)
    ffry, ffiy = fft.fft(Y)
    ffry, ffiy = ffry / len(Y), ffiy / len(Y)
    ffry, ffiy = ffry / (frq ** 2), ffiy / (frq ** 2)
    # Calculate the Quadrature
    # qst+i*qsi = ( xr - i*xi ) * (yr + i*yr)
    #           =  ( xr*yr + xi*yi) + i (-xi*yr + xr*yi)
    if flag == 1:
        A = 2 * (ffrx * ffry + ffix * ffiy)
    elif flag == 2:
        A = 2 * (ffrx * ffiy - ffix * ffry)
    A[int(len(A) / 2) :] = 0
    return A


def cal_blk(Ncut, frq, Pxx, Pyy, Pzz, Qxz, Qyz, Cxy):
    a1, b1 = zeros(len(Pxx)), zeros(len(Pxx))
    Mwd, Mds = zeros(len(Pxx)), zeros(len(Pxx))
    Hs = 4 * (num.sum(Pzz[Ncut:]) ** 0.5)
    Hsx = 4 * (num.sum(Pxx[Ncut:]) ** 0.5)
    Hsy = 4 * (num.sum(Pyy[Ncut:]) ** 0.5)
    Fp = frq[Ncut - 1 + num.argmax(Pzz[Ncut:])] / 6.28
    Fpx = frq[Ncut - 1 + num.argmax(Pxx[Ncut:])] / 6.28
    Fpy = frq[Ncut - 1 + num.argmax(Pyy[Ncut:])] / 6.28
    a1 = Qxz / (((Pxx + Pyy) * Pzz) ** 0.5)
    b1 = Qyz / (((Pxx + Pyy) * Pzz) ** 0.5)
    # a2 = (Pxx - Pyy) / (Pxx + Pyy)
    # b2 = (2 * Cxy) / (Pxx + Pyy)
    Mwd = vector.atan(b1 / a1)
    Mds = (2 * (1 - (a1 ** 2 + b1 ** 2) ** 0.5)) ** 0.5
    return Hs, Hsx, Hsy, Fp, Fpx, Fpy, Mwd, Mds


def iri(rb, Hs, Fp, Pzz, lat, lon, npx, Fil_Log) :
    rb.switch.value = True
    npx[0] = [npx.v, 0, 0]
    Nr = 1
    time.sleep(10)

    rb_sent = struct.pack("f", Hs)
    rb_sent += struct.pack("f", Fp)
    for j in range(len(Pzz)/4) :
        rb_sent += struct.pack("f", Pzz[j])
    rb_sent += struct.pack("f", lat)
    rb_sent += struct.pack("f", lon)
    rb.data_out = rb_sent
    status = rb.satellite_transfer()

    time.sleep(10)
    while status[0] > 8:
        Nr += 1
        rb.data_out = rb_sent   # Necessary ???
        status = rb.satellite_transfer()
        npx = LED(1, npx, [npx.v, 0, 0], [npx.v/3, npx.v/3, npx.v/3], 10)
    rb.switch.value = False
    return rb_sent, Nr
