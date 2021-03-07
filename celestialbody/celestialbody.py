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

def days_to_centuries(days):
    return days / 36525.


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
    # helps to read parameters file
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
        if name == "Earth":
            data_base_name = "EM Bary"
        else:
            data_base_name = name
        if category is None:
            if data_base_name in PLANETS:
                self._params = update_planet_params(data_base_name)
            elif data_base_name in ASTEROIDS:
                self._params = update_asteroid_params(data_base_name)
            elif data_base_name in COMETS:
                self._params = update_comet_params(data_base_name)
            else:
                print("Unknown object:", data_base_name)
                print("To check available bodies, see 'BODIES' list")
                return
        else:
            if category == "planet":
                self._params = update_planet_params(data_base_name)
                self.fullname = self._params["name"]
            elif category == "asteroid":
                self._params = update_asteroid_params(data_base_name)
                self.fullname = "("+str(self._params["number"])+") " + self._params["name"]
            elif category == "comet":
                self._params = update_comet_params(data_base_name)
                self.fullname = self._params["number"] + "/" + self._params["name"]
            else:
                print("Unknown category:", data_base_name)
                print("Available categories: 'planet', 'asteroid' or 'comet'")
                return

        self.tolerance = 1e-6 # tolerance for Kepler's equation numerical solver
        self.tolerance_rad = deg_to_rad(self.tolerance)
        self.max_iter = 1000

        self.name     = name
        self.category = self._params["category"]
        self.epoch    = self._params["epoch"]

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

        self._date = self.epoch
        self._a = self.a0
        self._e = self.e0
        self._i = self.i0
        self._L = self.L0
        self._varpi = self.varpi0
        self._Omega = self.Omega0

        self.n_pts_orbit = 1000  # nb of points for orbit list coordinates
        return

    @property
    def date(self):
        return self._date

    @date.setter
    def date(self, date):
        self._date = date

    @property
    def julian_ephemeris_date(self):
        return datetime_to_julian_ephemeris_date(self.date)

    @property
    def _days_since_epoch(self):
        t = self.julian_ephemeris_date
        t0 = datetime_to_julian_ephemeris_date(self.epoch)
        return t-t0

    @property
    def _centuries_since_j2000(self):
        return days_to_centuries(self.julian_ephemeris_date - j2000_jed)

    @property
    def period(self):
        """
        calulate period
        :return: period in days
        """
        return (2 * np.pi * np.sqrt((self.a * const.au) ** 3 / const.G / const.M_sun)).value / 24 / 3600

    @property
    def mean_motion(self):
        """
        mean motion in degree per day
        """
        return 360 / self.period

    @property
    def perihelion_passage_date(self):
        return self.date - datetime.timedelta(days=self.M / self.mean_motion)

    @property
    def a(self):
        if self.category == "planet":
            self._a = self.a0 + self.a_dot * self._centuries_since_j2000
        return self._a

    @property
    def e(self):
        if self.category == "planet":
            self._e = self.e0 + self.e_dot * self._centuries_since_j2000
        return self._e

    @property
    def i(self):
        if self.category == "planet":
            self._i = self.i0 + self.i_dot * self._centuries_since_j2000
        return self._i

    @property
    def L(self):
        if self.category == "planet":
            self._L = self.L0 + self.L_dot * self._centuries_since_j2000
        else:
            self._L = self.L0 + self.mean_motion * self._days_since_epoch
        return self._L

    @property
    def varpi(self):
        if self.category == "planet":
            self._varpi = self.varpi0 + self.varpi_dot * self._centuries_since_j2000
        return self._varpi

    @property
    def Omega(self):
        if self.category == "planet":
            self._Omega = self.Omega0 + self.Omega_dot * self._centuries_since_j2000
        return self._Omega

    @property
    def omega(self):
        """
        argument of perihelion in deg
        """
        return self.varpi-self.Omega

    @property
    def M(self):
        val = self.L - self.varpi
        if self.category == "planet":
            T = self._centuries_since_j2000
            corr = self.b * T**2 + self.c * np.cos(self.f*T) + self.s * np.sin(self.f*T)
            val += corr
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
                if np.abs(delta_E_rad) < self.tolerance_rad:
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

    @property
    def total_energy_per_kilogram(self):
        """
        :return: total energy per kilogram of the body
        """
        return (- const.G * const.M_sun / 2 / self.a).value

    @property
    def kinetic_energy_per_kilogram(self):
        """
        :return: kinetic energy per kilogram of the body
        """
        v = (np.sqrt(const.G * const.M_sun * (2 / self.r - 1 / self.a))).value
        return v ** 2 / 2

    @property
    def potential_energy_per_kilogram(self):
        """
        :return: total energy per kilogram of the body
        """
        return (- const.G * const.M_sun / self.r).value

    @property
    def area_constant(self):
        date = self.date
        self.date = self.perihelion_passage_date
        C = np.sqrt(2 * self.r**2 * self.kinetic_energy_per_kilogram)
        self.date = date
        return C

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

    @property
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

    @property
    def x(self):
        """
        :return: x coordinate in orbital plane at date
        """
        val, _, _ = self.orbital_heliocentric_coordinates
        return val

    @property
    def y(self):
        """
        :return: y coordinate in orbital plane at date
        """
        _, val, _ = self.orbital_heliocentric_coordinates
        return val

    @property
    def z(self):
        """
        :return: y coordinate in orbital plane at date
        """
        return 0

    @property
    def ecliptic_heliocentric_coordinates(self):
        """
        calculate ecliptic heliocentric coordinates at date
        :return: ecliptic heliocentric coordinates
        """
        x, y, z = self.orbital_heliocentric_coordinates
        return self.orbital_to_ecliptic_coordinates(x, y, z)

    @property
    def position(self):
        """
        just a simpler call to ecliptic_heliocentric_coordinates
        """
        return self.ecliptic_heliocentric_coordinates

    @property
    def X(self):
        """
        :return: x coordinate as ecliptic heliocentric coordinate at date
        """
        val, _, _ = self.position
        return val

    @property
    def Y(self):
        """
        :return: y coordinate as ecliptic heliocentric coordinate at date
        """
        _, val, _ = self.position
        return val

    @property
    def Z(self):
        """
        :return: z coordinate as ecliptic heliocentric coordinate at date
        """
        _, _, val = self.position
        return val

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

    def data(self, keyword, start=None, stop=None, step=None):
        date = self.date
        if start is None:
            start = datetime.datetime.today()
        if stop is None:
            stop = start + datetime.timedelta(days=self.period)
        if step is None:
            step = int((stop-start).days / 25)
        self.date = start
        tab, n = [], 0
        while self.date <= stop:
            if keyword == "days":
                val = n
            else:
                try:
                    prop = getattr(CelestialBody, keyword)
                    val = prop.fget(self)
                except:
                    print("Unknown property: "+keyword)
            tab.append(val)
            n += step
            self.date += datetime.timedelta(days=step)
        self.date = date
        return np.array(tab)

    def trajectory(self, start=None, stop=None, step=None):
        """
        :params:    start: datetime starting date, default is today
                    stop: datetime stop date, default is today + orbital period
                    step: time in days between two positions, default gives 25 points coordinates
        :return: x,y,z arrays of ecliptic heliocentric coordinates and n: array of days since start
        """
        positions = self.data("position", start=start, stop=stop, step=step)
        days = self.data("days", start=start, stop=stop, step=step)
        X,Y,Z = positions[:,0], positions[:,1], positions[:,2]
        return np.array(X), np.array(Y), np.array(Z), np.array(days)

    def data_position_txt(self, start=None, stop=None, step=None, filename=None, header=None, cols="xy", precision=5):
        """
        :params:    start: datetime starting date, default is today
                    stop: datetime stop date, default is today + orbital period
                    step: time in days between two positions, default gives 25 points coordinates
                    filename: name of the output file with .txt type, default is based on body name
                    header: header of the output file
                    cols: data shown in the output file, string with x, y, z and/or n in any order
                    precision: number of "significant digits", default is 5
        :return: filename
        """
        if start is None:
            start = datetime.datetime.today()
        if filename is None:
            filename = self.name.replace(" ", "_").lower() + ".txt"
        X, Y, Z, n = self.trajectory(start=start, stop=stop, step=step)
        if header is None:
            header = [u"# fichier " + filename + "\n",
                      u"####################################\n",
                      u"# days between two positions: %d \n" % (n[1]-n[0]),
                      u"# first line date: %2d/%2d/%4d\n" % (start.day, start.month, start.year),
                      u"####################################\n"
                     ]
        col_format = "{:<20}"
        col_label, col_data = "", []
        for col in cols:
            if col == "x":
                col_label += col_format.format("X (au)")
                col_data.append(X)
            elif col == "y":
                col_label += col_format.format("Y (au)")
                col_data.append(Y)
            elif col == "z":
                col_label += col_format.format("Z (au)")
                col_data.append(Z)
            elif col == "n":
                col_label += col_format.format("n (days)")
                col_data.append(n)

        with open(filename, "w") as f:
            f.writelines(header)
            f.write(col_label+"\n")
            for i in range(len(X)):
                line = ""
                for tab in col_data:
                    line += col_format.format(round(tab[i], precision-1))
                f.write(line+"\n")
        return filename