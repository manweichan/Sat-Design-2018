# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 16:27:00 2018

@author: wcgru
"""

from classDefinitions import Servicer, Thruster
import numpy as np
import matplotlib.pyplot as plt

# plt.ion()  # To immediately show plots

from astropy import units as u

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver


max_servicer_mass = 11000 # kg

dry_mass = 5000 # kg
wet_mass = max_servicer_mass - dry_mass 

numRatios = 50
ratioOfContinuousMassFuel = np.linspace(0.01,0.99,numRatios)

# set up parameters relevant to servicing mission profile
delI = 15# degrees 
# (total delta V for an average servicing mission)
serviceMissionDelV = 70 # m/s 
numberTrackingArray = np.zeros((numRatios,2))

# Set up orbit
initialOrbit = Orbit.circular(Earth, alt = 700e3 * u.m, inc = 98.6 * u.deg)

# Set up thrusters to use in the servicer
# based on  MONARC-445
Isp_imp = 235
thrust_imp  = 445
impulsiveThruster1 = Thruster(Isp = Isp_imp  , thrust = thrust_imp , engineMass = 1.6, power= 10)

impStr = 'Impulsive Thruster, Isp = ' +str (Isp_imp) + 's, thrust = ' + str(thrust_imp) + 'N'
mass_per_kw = 7 
# based on XIPS-13 Ion Engine
#==============================================================================
# xips_power = 0.42 #kw
# Isp_cont = 2500
# thrust_cont = 17.2 * 0.001
#==============================================================================

# based on XIPS-25 Ion Engine
xips_power = 2 #kw
Isp_cont = 3500
thrust_cont = 80* 0.001
contStr = 'Continuous Thruster, Isp = ' +str (Isp_cont) + 's, thrust = ' + str(thrust_cont) + 'N'
continuousThruster1 = Thruster(Isp = Isp_cont , thrust = thrust_cont, engineMass =mass_per_kw*xips_power , power= xips_power*1000)

# First analysis will be to compare the total number of missions and number of inclination changes based on fuel ratios
i = 0
for c_to_total_ratio in ratioOfContinuousMassFuel :
    impulsiveFuelMass = wet_mass * (1-c_to_total_ratio)
    continuousFuelMass = wet_mass * c_to_total_ratio
    servicer1 =  Servicer(initialOrbit, dryMass = dry_mass, \
                          impulsiveThruster = impulsiveThruster1, impulsiveFuelMass = impulsiveFuelMass , \
                          continuousThruster = continuousThruster1, continuousFuelMass = continuousFuelMass)
    
    numServiceMissions = 0
    numInclinationChanges = 0
    fuelRemaining = True
    continuousFuel = True
    impulsiveFuel = True
    counter = 0
    
    while fuelRemaining:
        if servicer1.execute_plane_change(delI,thrusterSelection = "continuous"):
            numInclinationChanges += 1
        else:
            continuousFuel = False

        if servicer1.execute_nominal_service_mission(serviceMissionDelV):
            numServiceMissions += 1
        else:
            impulsiveFuel  = False
            
        if not continuousFuel and not impulsiveFuel:
            fuelRemaining = False
            
        counter += 1
        #print('Number of fuel loops = ' + str(counter))

    numberTrackingArray[i,0] = numServiceMissions
    numberTrackingArray[i,1] = numInclinationChanges
    #print('On ratio loop: ' + str(i))
    i += 1

#%% Plots for Number of Missions Analysis

plt.plot(ratioOfContinuousMassFuel,numberTrackingArray)
plt.xlabel('Ratio of Continuous Thruster Fuel Mass to Total Fuel Mass')
plt.ylabel('Number of Missions')
plt.legend(['number of service missions ( ' + str(serviceMissionDelV) +' m/s each)','number of ' + str(delI) +'  deg inclination changes'])
plt.title('Starting Servicer Orbit: ' + str(servicer1.orbit) +'\n Total Servicer Mass: ' + str(max_servicer_mass)+' Total Fuel Mass: ' +str(wet_mass) + 'kg \n' + impStr + '\n' + contStr)

#%% More Detailed Analysis for Plane Changes
# next analysis will be to compare the amount of time / fuel / power it takes to execute a plane change for different continuous thrusters
delI = 15

# let's just take a 50% ratio between impulsive and continuous thrust fuel
# let's try a hall effect thruster (TsNIIMASH)

xips_power = 3.0 #kw
Isp_cont = 2000
thrust_cont = 1777* 0.001

continuousThruster1 = Thruster(Isp = Isp_cont , thrust = thrust_cont, engineMass =mass_per_kw*xips_power , power= xips_power*1000)

c_to_total_ratio = 0.5
impulsiveFuelMass = wet_mass * (1-c_to_total_ratio)
continuousFuelMass = wet_mass * c_to_total_ratio 
servicer1 =  Servicer(initialOrbit, dryMass = dry_mass, \
                          impulsiveThruster = impulsiveThruster1, impulsiveFuelMass = impulsiveFuelMass , \
                          continuousThruster = continuousThruster1, continuousFuelMass = continuousFuelMass)


(delV, time) = servicer1.execute_plane_change(delI,thrusterSelection = "continuous")

# Try using half of all impulsive thruster fuel first:
print('Total Delta-V required: ' + str(delV *0.001) + 'km/s')
print('Inclination Change of ' + str(delI) + ' deg, takes ' + str(time/(3600*24)) + ' days')


#%% Power is super important!
# make a plot of OAP vs. inclination change costs and time ( look for inflection points) 
# look at hall effect thruster (SME:SMAD), NEXT ion thruster (NASA paper), and other ion thrusters (SME:SMAD)


#%% Mission Based plots
# take a certain mission profile (X number refuels, Y numberr replacements, etc) and calculate revenue over time.

# should I make a separate mission class? or just have a servicer take in a mission one by one and have a mission directionary
