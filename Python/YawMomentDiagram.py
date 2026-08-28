import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as curve

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