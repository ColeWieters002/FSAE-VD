import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as curve
from Python.YawMomentDiagram.YMDSim import Tire, solve
from Python.Parameters import VehicleParameters as vp
from Parameters.VehicleParameters import LBF2N,N2LBF,FTLB2NM,NM2FTLB,FT2M,M2FT,IN2M,M2IN,RAD2DEG,DEG2RAD

#Inputs
Vx = 11.75 #m/s
beta = 0.0 * DEG2RAD
delta = 0.0 * DEG2RAD

r, Ay, Mz = solve(Vx, beta, delta, vp, tire)

print("r =", r)
print("Ay =", Ay)
print("Mz =", Mz)



# Sweep steering angle delta

    # Sweep sideslip beta

        # Define yaw rate r / curvature condition

        # For each wheel

            # Find wheel states
            #   wheel velocity components
            #   steer angle
            #   camber, etc.

            # Find vertical load Fz
            #   static load
            #   lateral load transfer
            #   aero if included

            # Find slip angle alpha

            # Find lateral force Fy

            # Find tire aligning moment Mz

        # Sum lateral force
        # Fy_total = Fy_FL + Fy_FR + Fy_RL + Fy_RR

        # Compute lateral acceleration
        # Ay = Fy_total / mass

        # Compute total vehicle yaw moment about CG
        # Mz_vehicle =
        #       moments caused by tire forces
        #     + tire aligning moments

        # Store:
        # Ay[i, j]
        # Mz[i, j]