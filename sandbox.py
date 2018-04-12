import numpy as np
import matplotlib.pyplot as plt

# plt.ion()  # To immediately show plots

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver


plt.style.use("seaborn")  # Recommended

ss_i = Orbit.circular(Earth, alt = 400e3 * u.m, inc = 30 * u.deg)

print(ss_i.v)

ss_im = ss_i.state.v.to(u.m/u.s)

print(ss_im)


r2d = 180/np.pi
d2r = np.pi/180

def delVIncChange(v,delI):
	#delI is in degrees
	delIrad = delI * d2r
	delV = 2*v*np.sin(delIrad/2)
	return delV

delVTest = delVIncChange(np.linalg.norm(ss_i.v),10)

print(np.linalg.norm(ss_i.v), 'km/s')

print(delVTest, 'km/s')


# print(ss_i)

# ss_f = ss_i.inc_change(inc = 40 * u.deg)
# ss_f = Maneuver.inclination(ss_i,inc = ss_i.inc + 10 * u.deg)


# print(ss_f.state)  # Updated after propagation
# print(ss_f.time_of_flight)  # Added up after propagation
# print(ss_f.total_cost)

