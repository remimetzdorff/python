import numpy as np
import datetime
import astropy.constants as const
import pandas
import os


j2000_date = datetime.datetime(2000, 1, 1, 12)
j2000_jed  = 2451545.


def update_planet_params(filename="p_elem_t2.txt"):
    """
    updates planet params from p_elem_t2.txt
    data are available at https://ssd.jpl.nasa.gov/txt/p_elem_t2.txt
    returns pandas Series with:
        index: planet names
        value: dict with keplerian parameters
    """
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
                      epoch                    = datetime.datetime(2000, 1, 1, 12),
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
        name = line_add_value[:10].strip(" ")
        if name in bodies.index:
            b = float(line_add_value[10:24].strip(" "))
            val = line_add_value[24:38].strip(" ")
            if len(val)>0:
                c = float(line_add_value[24:38].strip(" "))
                s = float(line_add_value[38:52].strip(" "))
                f = float(line_add_value[52:].strip(" "))
        planet = bodies[name]
        planet.update({"b":b, "c":c, "s":s, "f":f})
        n+=1
    return bodies


class CelestialBody():

    def __init__(self, name, verbose=False):
        self.verbose = verbose

        self._planets = update_planet_params()
        self._bodies = self._planets
        #self._bodies.append(self._comets)

        self.tolerance = 1e-6
        self.max_iter = 1000

        if name not in self._bodies.index:
            print("Unknown object:", name)
            print("Available objects:", self._bodies.index.values)
            return

        self._params        = self._bodies[name]
        self.category       = self._params["category"]
        self.epoch          = self._params["epoch"]
        self.name           = self._params["name"]
        if self.category == "planet":
            # Values
            self.a0     = self._params["semimajor_axis"]
            self.e0     = self._params["eccentricity"]
            self.i0     = self._params["inclination"]
            self.L0     = self._params["mean_longitude"]
            self.varpi0 = self._params["perihelion_longitude"]
            self.Omega0 = self._params["ascending_node_longitude"]
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

            self.date = self.epoch
            self.a = self.a0
            self.e = self.e0
            self.i = self.i0
            self.L = self.L0
            self.varpi = self.varpi0
            self.Omega = self.Omega0

            self.n_pts_orbit = 1000 # nb of points for orbit list coordinates
        else:
            print("Unknown object category:", self.category)
        return

    @property
    def date(self):
        return self.__date

    @date.setter
    def date(self, date):
        self.__date = date
        self.update_kepler_params()

    def date_to_julian_ephemeris_date(self, date):
        """
        converts datetime date to julian ephemeris date in days
        """
        delta = date - j2000_date
        return j2000_jed + delta.days + delta.seconds/24./3600.

    def days_to_centuries(self, n_days):
        return n_days / 36525.

    def update_kepler_params(self):
        if self.verbose:
            print("Updating, Kepler parameters")
        T_eph = self.date_to_julian_ephemeris_date(self.date)
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
            T_eph = self.date_to_julian_ephemeris_date(self.date)
            T = self.days_to_centuries(T_eph)
            corr = self.b * T**2 + self.c * np.cos(self.f*T) + self.s * np.sin(self.f*T)
        else:
            corr = 0
        val = self.L - self.varpi + corr
        val += 180
        return val % 360 - 180

    @property
    def E(self):
        """
        numerically solve M = E - e_star*sin(E)
        :return: solution E
        """
        val_M = self.M
        val_M_rad = val_M * np.pi / 180
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
        val_E = val_E_rad*180/np.pi
        return val_E

    def orbital_to_ecliptic_coordinates(self, x, y, z):
        omega_rad = self.omega * np.pi / 180
        Omega_rad = self.Omega * np.pi / 180
        i_rad     = self.i * np.pi / 180
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
        val_E_rad = self.E * np.pi / 180
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
        :return: x,y,z lists of ecliptic heliocentric coordinates
        """
        x, y, z = self.orbital_heliocentric_orbit()
        return self.orbital_to_ecliptic_coordinates(x, y, z)

    @property
    def orbit(self):
        """
        just a simpler call to ecliptic_heliocentric_orbit()
        """
        return self.ecliptic_heliocentric_orbit()