import os
import numpy as np
from scipy.interpolate import interp1d as curve

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_TIRE_DIR = os.path.join(_THIS_DIR, "Tires")

'''
UNIT CONVERSIONS
'''
LBF2N   = 4.4482216
N2LBF   = 1.0 / LBF2N
FTLB2NM = 1.3558
NM2FTLB = 1 / FTLB2NM
FT2M    = 0.3048
M2FT    = 1 / FT2M
IN2M    = 0.0254
M2IN    = 1 / IN2M
RAD2DEG = 180.0 / np.pi
DEG2RAD = np.pi / 180.0


'''
ACTUATION
'''
FrontSpringRate = 820 * LBF2N / IN2M #N/m
RearSpringRate = 640 * LBF2N / IN2M #N/m
#fmrWheelposistions = []
#fmrspringposistions = []
#FMRCurve = curve(fmrWheelposistions, fmrspringposistions, kind='quadratic')
FrontRollStiffness = 1000 #N*m/rad
RearRollStiffness = 1000 #N*m/rad

#LOOK AT PHONE FOR NUMEBRS

FrontRollStiffness = 10000 #N/m
RearRollStiffness = 10000 #N/m
#FrontHeaveStiffness
#RearHeaveStiffness
#ChassisStiffness
#LLTD
'''
#BRAKES
'''


'''
#HARDPOINTS
'''
##Mass/Inertia
Mass_kg = 176.0 # ~390 lbs
Driver_kg = 76.0
TotalMass_kg = Mass_kg + Driver_kg
UnsprungMass_kg = 42.0
SprungMass_kg = Mass_kg - UnsprungMass_kg # Assumed 42 kg unsprung
WeightDist = 0.5 # front
CG_mm = 324.0 # ~12.75 in
YawInertia = 92.0 # kg*m^2
# SprungRollInertia
# SprungPitchInertia

##WheelSpacing
Wheelbase_mm = 1535.0 # ~60.5 in
FTrackwidth_mm = 1200.0 # ~47 in
RTrackwidth_mm = 1168.0 # ~46 in

##Camber
CamberBounds = [0, -2, -3] # deg

_CamberCurve = curve([-25.4, 0, 25.4],CamberBounds,kind='quadratic',fill_value='extrapolate')

def Camber_By_Travel_deg(travel_mm, side):
    camber = float(_CamberCurve(travel_mm))

    if side.lower() == "left":
        return camber
    elif side.lower() == "right":
        return -camber
    else:
        raise ValueError("side must be 'left' or 'right'")

##RollCenter
FrontRollCenter_mm = -0.02 * IN2M
RearRollCenter_mm = 2.474 * IN2M 

##Steering
Ackerman = 0.0
FrontStaticToe_deg = 0.0
RearStaticToe_deg = 0.0

##Caster
FrontCaster_deg = 8.23
RearCaster_deg = 2.81

##KPI
FrontKPI_deg = 0.0
RearKPI_deg = 14.9

'''
#HUBS
'''


'''
#UPWRITES
'''



Gravity = 9.8 #m/s^2
AirDensity = 1.225 #kg/m^3
A = 1 #m^2
CL = -3.75
CD = 1.4
AeroBalance = .4 #front

#Will probably convert tire to class

TirePressure_bar = .827 #12psi
RollResistance = .4

#Tire Models
Hoosier_16x75_10_R20 = os.path.join(_TIRE_DIR, "Hoosier_16x75_10_R20.tir")
Hoosier_18x75_10_R20 = os.path.join(_TIRE_DIR, "Hoosier_18x75_10_R20.tir")
Hoosier_16x75_10_LC0 = os.path.join(_TIRE_DIR, "Hoosier_16x75_10_LC0.tir")
TireModel = Hoosier_16x75_10_R20
