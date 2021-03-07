CelestialBody is a project developed to calculate the position of an object in the solar system based on its keplerian parameters available online, for education purpose.

For planets, a guide is available to calculate heliocentric ecliptic J2000 coordinates : https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
More on heliocentric ecliptic J2000 coordinate : https://omniweb.gsfc.nasa.gov/coho/helios/coor_des.html

For comets and asteroids, see course at https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Book%3A_Celestial_Mechanics_(Tatum)/10%3A_Computation_of_an_Ephemeris/10.07%3A_Calculating_the_Position_of_a_Comet_or_Asteroid

For natural satellites, one might use orbital parameters from https://ssd.jpl.nasa.gov/?sat_elem to compute the rough estimate of any known satellite position.

Some data calculated from CelestialBody have been checked "by eye" using a solar system simulator: https://theskylive.com/3dsolarsystem
However, precise ephemeris computation is not guaranteed with this project. If high accuracy is required, you might need to check for data from the JPL HORIZONS system: https://ssd.jpl.nasa.gov/?horizons