#Vehicle Parameters
from scipy.interpolate import interp1d as curve

Mass_kg = 176.0 #~390lbs
CG_mm = 324.0 #~12.75in
Driver_kg = 76.0
Wheelbase_mm = 1535.0 #~60.5in
FTrackwidth_mm = 1200.0 #~47in
RTrackwidth_mm = 1168.0 #~46in
WeightDist = .5 #front
YawInertia = 92.0 #kg*m^2
SprungMass_kg = Mass_kg-42.0 #Assumed 42kg Unsprung
#SprungRollInertia
#SprungPitchInertia

FrontRollStiffness = 3000 #N/m
RearRollStiffness = 1000 #N/m
#FrontHeaveStiffness
#RearHeaveStiffness
#ChassisStiffness
#LLTD

#KPI
#CasterAngle

FrontStaticToe = 0 #deg
RearStaticToe = 0 #deg

#Front and Rear Camber will be separated later
CamberBounds = [0, -2, -3] #deg
Camber_By_Travel_deg = curve([-25.4, 0, 25.4], CamberBounds, kind='quadratic', fill_value='extrapolate') #mm

FrontRollCenter_mm = 0
RearRollCenter_mm = 50.8 #2in

Ackerman = 0.0

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
Hoosier_16x75_10_R20 = "Python\\Tires\\Hoosier_16x75_10_R20.tir"
Hoosier_18x75_10_R20 = "Hoosier_18x75_10_R20.tir"
Hoosier_16x75_10_LC0 = "Hoosier_16x75_10_LC0.tir"
TireModel = Hoosier_16x75_10_R20




