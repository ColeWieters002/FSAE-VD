import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import VehicleParameters as vp
from BicycleSim import Tire, solve, DEG2RAD
import time


def main():
    start = time.perf_counter()
    tire = Tire(vp.TireModel, vp.TirePressure_bar)
    print(f"Loaded tire: {tire.path}")

    #run that john
    Vx    = 10 #m/s ~25mph
    delta = 10 * DEG2RAD #rad

    result = solve(Vx, delta, vp, tire)
    beta, r, Ay = result
    end = time.perf_counter()
    Runtime = end-start
    print(f"Sideslip Angle = {beta * 57.296:.3f}, Yaw Rate = {r * 57.296:.3f}, Lateral Gs = {Ay / 9.8:.3f}")
    radius = (Vx**2)/Ay #m
    SkidpadTime = (2*np.pi)/r #sec
    print(f"radius= {radius:.3f}m")
    print(f"SkidpadTime= {SkidpadTime:.3f}s")
    print(f"Runtime = {Runtime:.3f}s")

if __name__ == "__main__":
    main()


