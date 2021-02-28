import numpy as np
import datetime
import astropy.constants as const
import pandas
import os


j2000_date = datetime.datetime(2000, 1, 1, 12)
j2000_jed  = 2451545.

planets_filename   = "p_elem_t2.txt"
asteroids_filename = "ELEMENTS.NUMBR"
comets_filename    = "ELEMENTS.COMET"

def rad_to_deg(angle):
    return angle * 180 / np.pi

def deg_to_rad(angle):
    return angle * np.pi / 180

def datetime_to_julian_ephemeris_date(date):
    """
    converts datetime date to julian ephemeris date in days
    """
    delta = date - j2000_date
    return j2000_jed + delta.days + delta.seconds / 24. / 3600.

def julian_ephemeris_date_to_datetime(jed):
    """
    converts datetime date to julian ephemeris date in days
    """
    delta = jed - j2000_jed
    return j2000_date + datetime.timedelta(days=int(delta),seconds=(delta-int(delta))*24*3600)

def update_planet_params(name):
    """
    updates planet params from p_elem_t2.txt
    data are available at https://ssd.jpl.nasa.gov/txt/p_elem_t2.txt
    returns pandas Series with:
        index: planet names
        value: dict with keplerian parameters
    """
    filename = planets_filename
    folder = os.path.dirname(__file__)+"/"
    with open(folder+filename) as file:
        content = file.readlines()

    # Main paramaters
    n_header, n_footer = 17, 18
    n = n_header
    names, planets = [], []
    while n < (len(content)-n_footer):
        line_value = content[n].strip("\n")
        line_rate  = content[n+1].strip("\n")
        planet = dict(
                      category                 = "planet",
                      epoch                    = j2000_date,
                      # Values
                      name                     = line_value[:9].strip(" "),
                      semimajor_axis           = float(line_value[9:20].strip(" ")),
                      eccentricity             = float(line_value[20:36].strip(" ")),
                      inclination              = float(line_value[36:52].strip(" ")),
                      mean_longitude           = float(line_value[52:70].strip(" ")),
                      perihelion_longitude     = float(line_value[70:86].strip(" ")),
                      ascending_node_longitude = float(line_value[86:].strip(" ")),
                      # Rates
                      semimajor_axis_rate           = float(line_rate[9:20].strip(" ")),
                      eccentricity_rate             = float(line_rate[20:36].strip(" ")),
                      inclination_rate              = float(line_rate[36:52].strip(" ")),
                      mean_longitude_rate           = float(line_rate[52:70].strip(" ")),
                      perihelion_longitude_rate     = float(line_rate[70:86].strip(" ")),
                      ascending_node_longitude_rate = float(line_rate[86:].strip(" "))
                     )
        # Additional parameters
        planet.update({"b":0, "c":0, "s":0, "f":0})
        names.append(planet["name"])
        planets.append(planet)
        n+=2
    bodies = pandas.Series(planets, index=names)

    # Update additional parameters for far planets
    n_header, n_footer = 47, 1
    n = n_header
    while n < (len(content)-n_footer):
        b, c, s, f = 0, 0, 0, 0
        line_add_value = content[n].strip("\n")
        current_name = line_add_value[:10].strip(" ")
        if current_name in bodies.index:
            b = float(line_add_value[10:24].strip(" "))
            val = line_add_value[24:38].strip(" ")
            if len(val)>0:
                c = float(line_add_value[24:38].strip(" "))
                s = float(line_add_value[38:52].strip(" "))
                f = float(line_add_value[52:].strip(" "))
        planet = bodies[current_name]
        planet.update({"b":b, "c":c, "s":s, "f":f})
        n+=1
    return bodies[name]

def planets_list():
    names = ["Mercury", "Venus", "EM Bary", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
    return names

def read_asteroid_params(row):
    """
        read a given row from "ELEMENTS.NUMBR"
    """
    num_asteroid = dict(
        category="asteroid",
        number=int(row[:6].strip(" ")),
        modified_epoch=int(row[25:30].strip(" ")),
        # Values
        name=row[7:24].strip(" "),
        semimajor_axis=float(row[31:41].strip(" ")),
        eccentricity=float(row[42:52].strip(" ")),
        inclination=float(row[53:62].strip(" ")),
        perihelion_argument=float(row[63:72].strip(" ")),
        ascending_node_longitude=float(row[73:82].strip(" ")),
        mean_anomaly=float(row[83:94].strip(" "))
    )
    return num_asteroid

def update_asteroid_params(name):
    """
        updates params from a specified numbered asteroid from ELEMENTS.NUMBR
        data are available at https://ssd.jpl.nasa.gov/?sb_elem
        return: dict with asteroid params
    """
    filename = asteroids_filename
    folder = os.path.dirname(__file__) + "/"
    with open(folder + filename) as file:
        content = file.readlines()
    for row in content:
        if name in row:
            break
    num_asteroid = read_asteroid_params(row)
    num_asteroid["epoch"] = julian_ephemeris_date_to_datetime(num_asteroid["modified_epoch"] + 2400000.5)
    num_asteroid["perihelion_longitude"] = num_asteroid["ascending_node_longitude"] + num_asteroid[
        "perihelion_argument"]
    num_asteroid["mean_longitude"] = num_asteroid["mean_anomaly"] + num_asteroid["perihelion_longitude"]
    return num_asteroid

def asteroids_list():
    filename = asteroids_filename
    folder = os.path.dirname(__file__) + "/"
    with open(folder + filename) as file:
        content = file.readlines()
    names = []
    for row in content[2:]:
        name = row[7:24].strip(" ")
        names.append(name)
    return names

def read_comets_params(row):
    """
        read a given row from "COMETS.NUMBR"
    """
    comet = dict(
        category="comet",
        number=row[:4].strip(" "),
        modified_epoch=int(row[44:51].strip(" ")),
        # Values
        name=row[5:43].strip(" "),
        perihelion_distance=float(row[52:63].strip(" ")),
        eccentricity=float(row[64:74].strip(" ")),
        inclination=float(row[75:84].strip(" ")),
        perihelion_argument=float(row[85:94].strip(" ")),
        ascending_node_longitude=float(row[95:104].strip(" ")),
        perihelion_passage_time=row[105:119].strip(" ")
    )
    return comet

def update_comet_params(name):
    """
        updates params from a specified numbered asteroid from COMETS.NUMBR
        data are available at https://ssd.jpl.nasa.gov/?sb_elem
        return: dict with comet params
    """
    filename = comets_filename
    folder = os.path.dirname(__file__) + "/"
    with open(folder + filename) as file:
        content = file.readlines()
    for row in content:
        if name in row:
            break
    comet = read_comets_params(row)
    comet["epoch"] = julian_ephemeris_date_to_datetime(comet["modified_epoch"] + 2400000.5)
    comet["perihelion_longitude"] = comet["ascending_node_longitude"] + comet["perihelion_argument"]
    comet["semimajor_axis"] = comet["perihelion_distance"] / (1 - comet["eccentricity"])
    date = comet["perihelion_passage_time"]
    year, month, day, frac_day = int(date[:4]), int(date[4:6]), int(date[6:8]), float(date[8:])
    comet["perihelion_passage_time"] = datetime.datetime(year, month, day) + datetime.timedelta(days=frac_day)
    t = datetime_to_julian_ephemeris_date(comet["epoch"])
    T = datetime_to_julian_ephemeris_date(comet["perihelion_passage_time"])
    P = (2 * np.pi * np.sqrt((comet["semimajor_axis"] * const.au) ** 3 / const.G / const.M_sun)).value / 24 / 3600
    comet["mean_anomaly"] = 360 / P * (t - T)
    comet["mean_longitude"] = comet["mean_anomaly"] + comet["perihelion_longitude"]
    return comet

def comets_list():
    filename = comets_filename
    folder = os.path.dirname(__file__) + "/"
    with open(folder + filename) as file:
        content = file.readlines()
    names = []
    for row in content[2:]:
        name = row[5:43].strip(" ")
        names.append(name)
    return names

def enumerate_string(content):
    print(content[0].strip("\n"))
    print(content[1].strip("\n"))
    line1, line2 = "", ""
    for i in range(len(content[0].strip("\n"))):
        line1+=(str(i%10))
        if i%10 == 0:
            dozens = str(i//10)
            while len(dozens)<10:
                dozens += " "
            line2 += dozens
    print(line1)
    print(line2[:len(line1)])
    return

PLANETS   = planets_list()
ASTEROIDS = asteroids_list()
COMETS    = comets_list()
BODIES    = PLANETS + ASTEROIDS + COMETS


class CelestialBody():

    def __init__(self, name, category=None, verbose=False):
        self.verbose = verbose

        if category is None:
            if name in PLANETS:
                self._params = update_planet_params(name)
            elif name in ASTEROIDS:
                self._params = update_asteroid_params(name)
            elif name in COMETS:
                self._params = update_comet_params(name)
            else:
                print("Unknown object:", name)
                print("To check available bodies, see 'BODIES' list")
                return
        else:
            if category == "planet":
                self._params = update_planet_params(name)
                self.fullname = self._params["name"]
            elif category == "asteroid":
                self._params = update_asteroid_params(name)
                self.fullname = "("+str(self._params["number"])+") " + self._params["name"]
            elif category == "comet":
                self._params = update_comet_params(name)
                self.fullname = self._params["number"] + "/" + self._params["name"]
            else:
                print("Unknown category:", name)
                print("Available categories: 'planet', 'asteroid' or 'comet'")
                return

        self.tolerance = 1e-6 # tolerance for Kepler's equation numerical solver
        self.max_iter = 1000

        self.name           = self._params["name"]
        self.category       = self._params["category"]
        self.epoch          = self._params["epoch"]

        self.a0 = self._params["semimajor_axis"]
        self.e0 = self._params["eccentricity"]
        self.i0 = self._params["inclination"]
        self.L0 = self._params["mean_longitude"]
        self.varpi0 = self._params["perihelion_longitude"]
        self.Omega0 = self._params["ascending_node_longitude"]

        if self.category == "planet":
            # Rates
            self.a_dot     = self._params["semimajor_axis_rate"]
            self.e_dot     = self._params["eccentricity_rate"]
            self.i_dot     = self._params["inclination_rate"]
            self.L_dot     = self._params["mean_longitude_rate"]
            self.varpi_dot = self._params["perihelion_longitude_rate"]
            self.Omega_dot = self._params["ascending_node_longitude_rate"]
            # Additional parameters
            self.b = self._params["b"]
            self.c = self._params["c"]
            self.s = self._params["s"]
            self.f = self._params["f"]

        self._old_date = self.epoch
        self._old_E = None

        self.date = self.epoch
        self.a = self.a0
        self.e = self.e0
        self.i = self.i0
        self.L = self.L0
        self.varpi = self.varpi0
        self.Omega = self.Omega0

        self.n_pts_orbit = 1000  # nb of points for orbit list coordinates
        return

    @property
    def date(self):
        return self.__date

    @date.setter
    def date(self, date):
        self.__date = date
        if self.category == "planet":
            self.update_kepler_params()

    def days_to_centuries(self, n_days):
        return n_days / 36525.

    def update_kepler_params(self):
        if self.verbose:
            print("Updating, Kepler parameters")
        T_eph = datetime_to_julian_ephemeris_date(self.date)
        T = self.days_to_centuries(T_eph-j2000_jed)
        self.a = self.a0 + self.a_dot * T
        self.e = self.e0 + self.e_dot * T
        self.i = self.i0 + self.i_dot * T
        self.L = self.L0 + self.L_dot * T
        self.varpi = self.varpi0 + self.varpi_dot * T
        self.Omega = self.Omega0 + self.Omega_dot * T
        return

    @property
    def period(self):
        """
        calulate period
        :return: period in days
        """
        P = (2 * np.pi * np.sqrt((self.a * const.au) ** 3 / const.G / const.M_sun)).value / 24 / 3600
        return P

    @property
    def omega(self):
        """
        argument of perihelion in deg
        """
        return self.varpi-self.Omega

    @property
    def M(self):
        if self.category == "planet":
            T_eph = datetime_to_julian_ephemeris_date(self.date)
            T = self.days_to_centuries(T_eph)
            corr = self.b * T**2 + self.c * np.cos(self.f*T) + self.s * np.sin(self.f*T)
            val = self.L - self.varpi + corr
        elif (self.category == "asteroid") or (self.category == "comet"):
            M0 = self.L0 - self.varpi0
            t = datetime_to_julian_ephemeris_date(self.date)
            t0 = datetime_to_julian_ephemeris_date(self.epoch)
            val = M0 + 360/self.period * (t-t0)
        val += 180
        return val % 360 - 180

    @property
    def E(self):
        """
        numerically solve Kepler's equation M = E - e_star*sin(E)
        :return: solution E
        """
        if (self.date != self._old_date) or (self._old_E is None):
            val_M = self.M
            val_M_rad = deg_to_rad(val_M)
            # Newton's method
            val_E_rad = val_M_rad + self.e * np.sin(val_M_rad)
            success = False
            for i in range(self.max_iter):
                delta_M_rad = val_M_rad - (val_E_rad - self.e * np.sin(val_E_rad))
                delta_E_rad = delta_M_rad / (1 - self.e * np.cos(val_E_rad))
                val_E_rad += delta_E_rad
                if np.abs(delta_E_rad) < self.tolerance*np.pi/180:
                    success = True
                    break
            if not success:
                print("FAIL TO REACH CLOSE ENOUGH VALUE")
            val_E = rad_to_deg(val_E_rad)
            self._old_E = val_E
            self._old_date = self.date
        else:
            val_E = self._old_E
        return val_E

    @property
    def nu(self):
        E_rad = deg_to_rad(self.E)
        val = np.sqrt((1+self.e)/(1-self.e)) * np.tan(E_rad/2)
        return rad_to_deg(2*np.arctan(val))

    @property
    def r(self):
        return self.a * (1 - self.e**2) / (1+ self.e * np.cos(deg_to_rad(self.nu)))

    def orbital_to_ecliptic_coordinates(self, x, y, z):
        omega_rad = deg_to_rad(self.omega)
        Omega_rad = deg_to_rad(self.Omega)
        i_rad     = deg_to_rad(self.i)
        X = (np.cos(omega_rad) * np.cos(Omega_rad) - np.sin(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * x \
            + (-np.sin(omega_rad) * np.cos(Omega_rad) - np.cos(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * y
        Y = (np.cos(omega_rad) * np.sin(Omega_rad) + np.sin(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * x \
            + (-np.sin(omega_rad) * np.sin(Omega_rad) + np.cos(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * y
        Z = (np.sin(omega_rad) * np.sin(i_rad)) * x + (np.cos(omega_rad) * np.sin(i_rad)) * y
        return X, Y, Z

    def orbital_heliocentric_coordinates(self):
        """
        calculate coordinates in orbital plane at date
        :return: x, y, z
        """
        val_E_rad = deg_to_rad(self.E)
        x = self.a * (np.cos(val_E_rad)-self.e)
        y = self.a * np.sqrt(1 - self.e**2) * np.sin(val_E_rad)
        z = 0
        return x, y, z

    def ecliptic_heliocentric_coordinates(self):
        """
        calculate ecliptic heliocentric coordinates at date
        :return: ecliptic heliocentric coordinates
        """
        x, y, z = self.orbital_heliocentric_coordinates()
        return self.orbital_to_ecliptic_coordinates(x, y, z)

    @property
    def position(self):
        """
        just a simpler call to ecliptic_heliocentric_coordinates()
        """
        return self.ecliptic_heliocentric_coordinates()

    def orbital_heliocentric_orbit(self):
        """
        gives orbit trace as x,y,z lists of orbital heliocentric coordinates along an entire revolution
        spacing between position is not related to actual speed of the body
        :return: x,y,z lists of orbital heliocentric coordinates
        """
        val_E_rad = np.linspace(-np.pi, np.pi, self.n_pts_orbit+1)
        x = self.a * (np.cos(val_E_rad) - self.e)
        y = self.a * np.sqrt(1 - self.e ** 2) * np.sin(val_E_rad)
        z = 0
        return x, y, z

    def ecliptic_heliocentric_orbit(self):
        """
        gives orbit trace as x,y,z lists of ecliptic heliocentric coordinates along an entire revolution
        :return: x,y,z arrays of ecliptic heliocentric coordinates
        """
        x, y, z = self.orbital_heliocentric_orbit()
        return self.orbital_to_ecliptic_coordinates(x, y, z)

    @property
    def orbit(self):
        """
        just a simpler call to ecliptic_heliocentric_orbit()
        """
        return self.ecliptic_heliocentric_orbit()

    def trajectory(self, start, stop, step):
        """
        :params:    start: datetime starting date
                    stop: datetime stop date
                    step: time in days between two positions
        :return: x,y,z arrays of ecliptic heliocentric coordinates and n: array of days since start
        """
        date = self.date
        self.date = start
        X,Y,Z,n = [], [], [], []
        while self.date <= stop:
            n.append((self.date-start).days)
            x,y,z = self.position
            X.append(x)
            Y.append(y)
            Z.append(z)
            self.date += datetime.timedelta(days=step)
        self.date = date
        return np.array(X), np.array(Y), np.array(Z), np.array(n)