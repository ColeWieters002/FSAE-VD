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

    # -----------------------------------------
    # Skidpad inputs
    # -----------------------------------------
    delta = 10 * DEG2RAD       # rad
    Vx_min = 1.0               # m/s
    Vx_max = 30.0              # m/s
    Vx_step = 0.1              # m/s

    # -----------------------------------------
    # Find highest valid speed
    # -----------------------------------------
    max_valid_Vx = None
    max_result = None

    Vx_values = np.arange(Vx_min, Vx_max + Vx_step, Vx_step)

    for Vx in Vx_values:

        try:
            beta, r, Ay = solve(Vx, delta, vp, tire)

        except RuntimeError as e:
            print(f"Solver failed at Vx = {Vx:.2f} m/s")
            print(e)
            break

        # -----------------------------------------
        # Check that the result is physically valid
        # -----------------------------------------

        if not np.isfinite(beta):
            print(f"Invalid beta at Vx = {Vx:.2f} m/s")
            break

        if not np.isfinite(r):
            print(f"Invalid yaw rate at Vx = {Vx:.2f} m/s")
            break

        if not np.isfinite(Ay):
            print(f"Invalid lateral acceleration at Vx = {Vx:.2f} m/s")
            break

        if Ay <= 0:
            print(f"Invalid lateral acceleration at Vx = {Vx:.2f} m/s")
            break

        # -----------------------------------------
        # Calculate radius and time
        # -----------------------------------------

        radius = Vx**2 / Ay

        if r <= 0:
            print(f"Invalid yaw rate at Vx = {Vx:.2f} m/s")
            break

        skidpad_time = (2 * np.pi) / r

        # -----------------------------------------
        # Save this as the highest valid solution
        # -----------------------------------------

        max_valid_Vx = Vx
        max_result = (beta, r, Ay, radius, skidpad_time)

        print(
            f"Vx = {Vx:6.2f} m/s | "
            f"Beta = {beta * 57.296:7.3f} deg | "
            f"r = {r * 57.296:8.3f} deg/s | "
            f"Ay = {Ay / 9.8:6.3f} G | "
            f"R = {radius:7.3f} m | "
            f"Time = {skidpad_time:7.3f} s"
        )

    # -----------------------------------------
    # Final result
    # -----------------------------------------

    if max_valid_Vx is not None:

        beta, r, Ay, radius, skidpad_time = max_result

        print("\n==============================")
        print("Highest Valid Solution")
        print("==============================")

        print(f"Velocity       = {max_valid_Vx:.3f} m/s")
        print(f"Velocity       = {max_valid_Vx * 2.23694:.3f} mph")
        print(f"Sideslip       = {beta * 57.296:.3f} deg")
        print(f"Yaw Rate       = {r * 57.296:.3f} deg/s")
        print(f"Lateral Gs     = {Ay / 9.8:.3f} G")
        print(f"Radius         = {radius:.3f} m")
        print(f"Skidpad Time   = {skidpad_time:.3f} s")

    else:
        print("\nNo valid solution found.")

    end = time.perf_counter()
    Runtime = end - start

    print(f"Runtime        = {Runtime:.3f} s")


if __name__ == "__main__":
    main()