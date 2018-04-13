from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver

from astropy import units as u

import numpy as np
import matplotlib.pyplot as plt


def inc_change_delV(orbit1, orbit2):
	inc1 = orbit1.inc
	inc2 = orbit2.inc
	delI = inc2 - inc1
	delIrad = delI * np.pi/180
	print(delIrad)
	v = np.linalg.norm(orbit1.v)
	delV = 2*v*np.sin(delIrad.value/2)
	return delV

#Test
ss_i = Orbit.circular(Earth, alt = 400 * u.km, inc = 30 * u.deg)
ss_i2 = Orbit.circular(Earth, alt = 400 * u.km, inc = 40 * u.deg)

print(inc_change_delV(ss_i,ss_i2),' km/s')