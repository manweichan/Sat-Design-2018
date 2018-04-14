from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver

from astropy import units as u

import numpy as np
import matplotlib.pyplot as plt



def delVHohmann(alt1, alt2):
	#alt1 is the altitude of the first orbit, and alt2 is the altitude
	#of the second orbit. 
	#Assumes circular orbit
	re = 6371 #radius of the earth in km
	ss_i = Orbit.circular(Earth, alt = alt1 * u.km)
	hoh = Maneuver.hohmann(ss_i, (re + alt2) * u.km)
	return hoh.get_total_cost(), hoh

#test
delV, hoh = delVHohmann(700, 36000-6371)
print('Delta V:', delV.value, 'k m/s')

print(hoh.impulses[0])
print(hoh.impulses[1])