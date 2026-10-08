import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import VehicleParameters as vp
from BicycleSim import Tire, solve, DEG2RAD
import time



def main():
    start=time.perf_counter()

    tire=Tire(vp.TireModel,vp.TirePressure_bar)
    print(f"Loaded tire: {tire.path}")

    Vx=11.75


    #Solve 6 degrees
    delta=17*DEG2RAD

    beta,r,Ay,phi,data=solve(
        Vx,
        delta,
        vp,
        tire,
        return_data=True
    )

    print("\n6 DEG")
    print(f"Sideslip Angle = {beta*57.296:.3f} deg")
    print(f"Yaw Rate = {r*57.296:.3f} deg/s")
    print(f"Lateral Gs = {Ay/9.8:.3f}")

    x0=np.array([beta,r,Ay,phi])


    #Use 6 degree solution as starting point for 7 degrees
    delta=30*DEG2RAD

    beta,r,Ay,phi,data=solve(
        Vx,
        delta,
        vp,
        tire,
        x0=x0,
        return_data=True
    )

    print("\n7 DEG")
    print(f"Sideslip Angle = {beta*57.296:.3f} deg")
    print(f"Yaw Rate = {r*57.296:.3f} deg/s")
    print(f"Lateral Gs = {Ay/9.8:.3f}")

    radius=Vx**2/abs(Ay)
    SkidpadTime=2*np.pi/abs(r)

    print(f"radius = {radius:.3f} m")
    print(f"SkidpadTime = {SkidpadTime:.3f} s")

    end=time.perf_counter()
    print(f"Runtime = {end-start:.3f} s")


if __name__=="__main__":
    main()