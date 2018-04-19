# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:22:46 2018

@author: wcgru
"""
# Servicer class defintion
#==============================================================================
# The goal is to build a class that is initialized with a certain:
#     dry mass
#     ImpulsiveThruster (another class that has a certain Isp, thrust, etc)
#     hydrazine fuel mass (for impulsive manuevers)
#     ContinuousThruster (another class that has certain Isp, thrust, etc)
#     low-thrust mass (for plane changes)
#     
#     then have methods that can use the impulsive thruster and continuous thruster for certain manuevers and 
#     then the fuel profile can be viewed over a "lifecycle"
#==============================================================================
    
import numpy as np
class Servicer:
    def __init__(self, orbit, dryMass, impulsiveThruster, impulsiveFuelMass, continuousThruster, continuousFuelMass):
        self.orbit = orbit
        self.dryMass = dryMass + impulsiveThruster.mass + continuousThruster.mass
        self.impulsiveThruster = impulsiveThruster
        self.impulsiveFuelMass = impulsiveFuelMass
        self.continuousThruster = continuousThruster
        self.continuousFuelMass = continuousFuelMass
        
        # use an orbit object from poliastro
        self.currentOrbit = orbit
        self.totalMass = self.dryMass + impulsiveFuelMass + continuousFuelMass
    
    def execute_nominal_service_mission(self, serviceMissionDelV):
        # modifies total mass and impulsive fuel mass of servicer object
        # returns True if mission is executed and False if not execute (not enough fuel)
        # DOES NOT modify orbit, assumes it leaves from and returns to the same orbit
        
        # use impulsive thrusters for this mission, calculate mass lost
        finalMass = self.impulsiveThruster.calc_final_mass(delV = serviceMissionDelV, initialMass = self.totalMass)
        
        # decrement impulsive fuel
        tempImpulsiveFuelMass = self.impulsiveFuelMass - (self.totalMass - finalMass)
        if tempImpulsiveFuelMass < 0:
            #print('Not enough fuel to complete mission')
            return False
        else:
            # update mass allocations
            self.impulsiveFuelMass = tempImpulsiveFuelMass
            self.totalMass = finalMass
            return True
    
    # need to capture fuel use from station-keeping and atttiude manuevers, maybe use a constant loss rate?
    
    def execute_plane_change(self,delI,thrusterSelection):
        # (future refactor?) takes in an orbit object (new_orbit) from poliAstro and calculates the delV required depending on the thruster selection
        V = np.linalg.norm(self.orbit.v) * 1000 # convert from km/s to m/s
        if thrusterSelection == "continuous":
            # using a continuous thruster uses Edelbaum’s analytic solution
            # beta_0 is the initial thrust vector yaw angle
            # same orbital velocity
            Vi = V
            Vf = V
            delI_rad = np.pi/180 * delI
            beta_0 = np.arctan2(np.sin(np.pi*0.5*delI_rad),Vi/Vf - np.cos(np.pi*0.5*delI_rad))
            delV_total = np.sqrt(Vi**2 - 2*Vi*Vf*np.cos(np.pi*0.5*delI_rad) + Vf**2)
            thrustAcceleration = self.continuousThruster.thrust/self.totalMass
            time_total = delV_total/thrustAcceleration
            
            finalMass = self.continuousThruster.calc_final_mass(delV = delV_total, initialMass = self.totalMass)
            
             # decrement continuous fuel
            tempContinuousFuelMass = self.continuousFuelMass - (self.totalMass - finalMass)
             
            if tempContinuousFuelMass < 0:
                #print('Not enough fuel to complete plane change')
                return False
            else:
                self.continuousFuelMass = tempContinuousFuelMass
                self.totalMass = finalMass
                return (delV_total, time_total)
             
        elif thrusterSelection == "impulsive":
            
            #delI is in degrees
            delI_rad = np.pi/180 * delI
            delV = 2*V*np.sin(delI_rad/2)
            finalMass = self.impulsiveThruster.calc_final_mass(delV = delV,  initialMass = self.totalMass)
        
            # decrement impulsive fuel
            tempImpulsiveFuelMass = self.impulsiveFuelMass - (self.totalMass - finalMass)
            if tempImpulsiveFuelMass < 0:
                #print('Not enough fuel to complete plane change')
                return False
            else:
                # update mass allocations
                self.impulsiveFuelMass = tempImpulsiveFuelMass
                self.totalMass = finalMass
                return (delV, 0)
        else:
            print('only continuous and impulsive options available')
            return None
        
    
# look at table 18-5 in SME:SMAD For hydrazine thrusters and tables 18-12/18-13 for continuous thrusters
class Thruster:
    def __init__(self, Isp, thrust, engineMass, power):
        self.Isp = Isp #(s)
        self.thrust = thrust # (N)
        self.mass = engineMass
        self.power = power
        
    def calc_final_mass(self, delV, initialMass):
        return initialMass/np.exp(delV/(self.Isp*9.82))