import numpy as np
import matplotlib.pyplot as plt

import VehicleParameters as vp
from VehicleParameters import DEG2RAD, NM2FTLB
from YMDSim import Tire, solve


# ==========================================
# SETTINGS
# ==========================================

tire = Tire(vp.TireModel,vp.TirePressure_bar)

# Speed sweep
Vx_values = np.arange(8, 20, 1.0)

# YMD sweep used to find max Ay at each speed
beta_values = np.arange(-12, 13, 1)
delta_values = np.arange(-12, 13, 1)


# ==========================================
# STORAGE
# ==========================================

stability_values = []
control_values = []
CSR_values = []
Ay_max_values = []
Mz_at_Ay_max_values = []


# ==========================================
# SPEED SWEEP
# ==========================================

for Vx in Vx_values:

    print(f"\nRunning Vx = {Vx:.1f} m/s")


    # --------------------------------------
    # CENTER DERIVATIVES
    # --------------------------------------

    beta_step = 1.0 * DEG2RAD
    delta_step = 1.0 * DEG2RAD


    # beta derivative
    Ay_bp, Mz_bp, phi_bp,_ = solve(
        Vx,
        +beta_step,
        0.0,
        vp,
        tire
    )

    Ay_bm, Mz_bm, phi_bm,_ = solve(
        Vx,
        -beta_step,
        0.0,
        vp,
        tire
    )

    dMz_dbeta = (
        Mz_bp - Mz_bm
    ) / (2.0 * beta_step)


    # delta derivative
    Ay_dp, Mz_dp, phi_dp,_ = solve(
        Vx,
        0.0,
        +delta_step,
        vp,
        tire
    )

    Ay_dm, Mz_dm, phi_dm,_  = solve(
        Vx,
        0.0,
        -delta_step,
        vp,
        tire
    )

    dMz_ddelta = (
        Mz_dp - Mz_dm
    ) / (2.0 * delta_step)


    # Convert from N*m/rad to ft-lb/rad
    dMz_dbeta_ftlb = dMz_dbeta * NM2FTLB
    dMz_ddelta_ftlb = dMz_ddelta * NM2FTLB


    # CSR
    if abs(dMz_dbeta) > 1e-9:
        CSR = dMz_ddelta / dMz_dbeta
    else:
        CSR = np.nan


    # Store
    stability_values.append(dMz_dbeta_ftlb)
    control_values.append(dMz_ddelta_ftlb)
    CSR_values.append(CSR)


    # --------------------------------------
    # FULL YMD FOR MAX AY
    # --------------------------------------

    Ay_max = -np.inf
    Mz_at_Ay_max = np.nan

    for delta_deg in delta_values:

        delta = delta_deg * DEG2RAD

        for beta_deg in beta_values:

            beta = beta_deg * DEG2RAD

            Ay, Mz, phi, _ = solve(
                Vx,
                beta,
                delta,
                vp,
                tire
            )

            Ay_g = Ay / vp.Gravity

            if Ay_g > Ay_max:

                Ay_max = Ay_g
                Mz_at_Ay_max = Mz


    Ay_max_values.append(Ay_max)

    Mz_at_Ay_max_values.append(
        Mz_at_Ay_max * NM2FTLB
    )


    print(
        f"Stability = {dMz_dbeta_ftlb:.2f} ft-lb/rad | "
        f"Control = {dMz_ddelta_ftlb:.2f} ft-lb/rad | "
        f"CSR = {CSR:.2f} | "
        f"Max Ay = {Ay_max:.2f} g"
    )


# Convert to numpy arrays
stability_values = np.array(stability_values)
control_values = np.array(control_values)
CSR_values = np.array(CSR_values)
Ay_max_values = np.array(Ay_max_values)
Mz_at_Ay_max_values = np.array(Mz_at_Ay_max_values)


# ==========================================
# PLOT 1 — STABILITY DERIVATIVE
# ==========================================

plt.figure(figsize=(10, 6))

plt.plot(
    Vx_values,
    stability_values,
    marker="o"
)

plt.axhline(0)

plt.xlabel("Vehicle Speed (m/s)")
plt.ylabel(r"$dM_z/d\beta$ (ft-lb/rad)")
plt.title("Center Stability Derivative vs Vehicle Speed")

plt.grid(True)
plt.tight_layout()


# ==========================================
# PLOT 2 — CONTROL DERIVATIVE
# ==========================================

plt.figure(figsize=(10, 6))

plt.plot(
    Vx_values,
    control_values,
    marker="o"
)

plt.xlabel("Vehicle Speed (m/s)")
plt.ylabel(r"$dM_z/d\delta$ (ft-lb/rad)")
plt.title("Center Control Derivative vs Vehicle Speed")

plt.grid(True)
plt.tight_layout()


# ==========================================
# PLOT 3 — CONTROL / STABILITY RATIO
# ==========================================

plt.figure(figsize=(10, 6))

plt.plot(
    Vx_values,
    CSR_values,
    marker="o"
)

plt.axhline(0)

plt.xlabel("Vehicle Speed (m/s)")
plt.ylabel("Control / Stability Ratio")
plt.title("Control-to-Stability Ratio vs Vehicle Speed")

plt.grid(True)
plt.tight_layout()


# ==========================================
# PLOT 4 — MAXIMUM LATERAL ACCELERATION
# ==========================================

plt.figure(figsize=(10, 6))

plt.plot(
    Vx_values,
    Ay_max_values,
    marker="o"
)

plt.xlabel("Vehicle Speed (m/s)")
plt.ylabel(r"Maximum $A_y$ (g)")
plt.title("Maximum Lateral Acceleration vs Vehicle Speed")

plt.grid(True)
plt.tight_layout()


# ==========================================
# OPTIONAL PLOT 5 — MZ AT MAX AY
# ==========================================

plt.figure(figsize=(10, 6))

plt.plot(
    Vx_values,
    Mz_at_Ay_max_values,
    marker="o"
)

plt.axhline(0)

plt.xlabel("Vehicle Speed (m/s)")
plt.ylabel(r"$M_z$ at Maximum $A_y$ (ft-lb)")
plt.title("Yaw Moment at Maximum Lateral Acceleration")

plt.grid(True)
plt.tight_layout()


plt.show()