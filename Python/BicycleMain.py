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
    Vx    = 45 * vp.MPH2MPS # m/s
    delta = 15.7 * DEG2RAD #rad

    beta, r, Ay, phi, data = solve(Vx, delta, vp, tire, return_data=True)
    end = time.perf_counter()
    Runtime = end-start
    print(f"Inputs: Vx = {Vx * vp.MPS2MPH:.3f}, delta = {delta / DEG2RAD}")
    print(f"Sideslip Angle = {beta * 57.296:.3f} (deg), Yaw Rate = {r * 57.296:.3f} (deg/s), Lateral Gs = {Ay / vp.Gravity:.3f}")
    print(f"Body Roll = {phi / DEG2RAD:.3f} (deg)")
    radius = (Vx**2)/Ay #m
    SkidpadTime = (2*np.pi)/r #sec
    print(f"radius= {radius:.3f}m")
    print(f"SkidpadTime= {SkidpadTime:.3f}s")
    print(f"Runtime = {Runtime:.3f}s")

    print("Camber (deg; physical wheel-local, negative=negative camber):")
    print(
        "  " + ", ".join(
            f"{corner}={data[f'Camber_{corner}_deg']:+.2f}"
            for corner in ("FL", "FR", "RL", "RR")
        )
    )
    print("Wheel travel (mm; +compression, -extension):")
    print(
        "  " + ", ".join(
            f"{corner}={data[f'Travel_{corner}_mm']:+.2f}"
            for corner in ("FL", "FR", "RL", "RR")
        )
    )

if __name__ == "__main__":
    main()


