import numpy as np
import matplotlib.pyplot as plt
import VehicleParameters as vp
from VehicleParameters import DEG2RAD
from YMDSim import Tire, solve

# ===================================
# SETUP
# ===================================

Vx = 11.75
tracks = np.arange(42, 52.1, 1)

beta_coarse = np.arange(-12, 12.1, 2)
delta_coarse = np.arange(-12, 12.1, 2)

FineStep = 0.25
FineWindow = 1.0

tire = Tire(vp.TireModel, vp.TirePressure_bar)

F0 = vp.FTrackwidth_mm
R0 = vp.RTrackwidth_mm

F0_in = F0 / 1000 / vp.IN2M
R0_in = R0 / 1000 / vp.IN2M

n = len(tracks)

AyMax = np.zeros((n,n))
CSR = np.zeros((n,n))
MzAy = np.zeros((n,n))
dBeta = np.zeros((n,n))
dDelta = np.zeros((n,n))


# ===================================
# SWEEP
# ===================================

for i, rear in enumerate(tracks):
    for j, front in enumerate(tracks):

        vp.FTrackwidth_mm = front * vp.IN2M * 1000
        vp.RTrackwidth_mm = rear * vp.IN2M * 1000

        # CENTER DERIVATIVES
        _, Mp, _, _ = solve(Vx,  DEG2RAD, 0, vp, tire, False)
        _, Mm, _, _ = solve(Vx, -DEG2RAD, 0, vp, tire, False)

        dB = (Mp-Mm)/(2*DEG2RAD)

        _, Mp, _, _ = solve(Vx, 0,  DEG2RAD, vp, tire, False)
        _, Mm, _, _ = solve(Vx, 0, -DEG2RAD, vp, tire, False)

        dD = (Mp-Mm)/(2*DEG2RAD)

        dBeta[i,j] = dB * vp.NM2FTLB
        dDelta[i,j] = dD * vp.NM2FTLB
        CSR[i,j] = dD/dB

        # COARSE AY SEARCH
        bestAy = -np.inf
        bestB = 0
        bestD = 0

        for d in delta_coarse:
            for b in beta_coarse:

                Ay, Mz, _, _ = solve(
                    Vx,
                    b*DEG2RAD,
                    d*DEG2RAD,
                    vp,
                    tire,
                    False
                )

                if Ay > bestAy:
                    bestAy = Ay
                    bestB = b
                    bestD = d

        # LOCAL FINE SEARCH
        beta_fine = np.arange(
            max(-12, bestB-FineWindow),
            min(12, bestB+FineWindow)+FineStep/2,
            FineStep
        )

        delta_fine = np.arange(
            max(-12, bestD-FineWindow),
            min(12, bestD+FineWindow)+FineStep/2,
            FineStep
        )

        bestAy = -np.inf
        bestMz = 0

        for d in delta_fine:
            for b in beta_fine:

                Ay, Mz, _, _ = solve(
                    Vx,
                    b*DEG2RAD,
                    d*DEG2RAD,
                    vp,
                    tire,
                    False
                )

                if Ay > bestAy:
                    bestAy = Ay
                    bestMz = Mz

        AyMax[i,j] = bestAy / vp.Gravity
        MzAy[i,j] = bestMz * vp.NM2FTLB

        print(
            f"F={front:.0f} R={rear:.0f} | "
            f"Ay={AyMax[i,j]:.4f}g | "
            f"CSR={CSR[i,j]:.2f}"
        )


# ===================================
# RESTORE
# ===================================

vp.FTrackwidth_mm = F0
vp.RTrackwidth_mm = R0

F,R = np.meshgrid(tracks,tracks)


# ===================================
# CONTOUR FUNCTION
# ===================================

def contour(Z, title, label):

    plt.figure(figsize=(9,7))

    c = plt.contourf(F,R,Z,20)
    plt.colorbar(c,label=label)

    plt.scatter(
        F0_in,
        R0_in,
        marker="x",
        s=100,
        label="Current"
    )

    plt.xlabel("Front Track Width (in)")
    plt.ylabel("Rear Track Width (in)")
    plt.title(title)

    plt.legend()
    plt.tight_layout()
    plt.show()


# ===================================
# CONTOURS
# ===================================

contour(
    AyMax,
    "Front / Rear Track vs Maximum Lateral Acceleration",
    "Maximum Ay (g)"
)

contour(
    CSR,
    "Front / Rear Track vs Control / Stability Ratio",
    "CSR"
)

contour(
    MzAy,
    "Front / Rear Track vs Mz at Maximum Ay",
    "Mz at Ay Max (ft-lb)"
)


# ===================================
# FRONT VS REAR EFFECT
# ===================================

iR = np.argmin(abs(tracks-R0_in))
jF = np.argmin(abs(tracks-F0_in))

plt.figure(figsize=(9,6))

plt.plot(
    tracks,
    AyMax[iR,:],
    "o-",
    label=f"Front Sweep (Rear={tracks[iR]:.0f} in)"
)

plt.plot(
    tracks,
    AyMax[:,jF],
    "s-",
    label=f"Rear Sweep (Front={tracks[jF]:.0f} in)"
)

plt.axvline(F0_in, linestyle="--")
plt.xlabel("Track Width (in)")
plt.ylabel("Maximum Ay (g)")
plt.title("Front vs Rear Track Width Effect on Ay")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# ===================================
# CSR FRONT VS REAR
# ===================================

plt.figure(figsize=(9,6))

plt.plot(
    tracks,
    CSR[iR,:],
    "o-",
    label="Front Track Sweep"
)

plt.plot(
    tracks,
    CSR[:,jF],
    "s-",
    label="Rear Track Sweep"
)

plt.xlabel("Track Width (in)")
plt.ylabel("Control / Stability Ratio")
plt.title("Front vs Rear Track Width Effect on CSR")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# ===================================
# DERIVATIVES AT CURRENT REAR TRACK
# ===================================

fig,ax1 = plt.subplots(figsize=(9,6))
ax2 = ax1.twinx()

l1 = ax1.plot(
    tracks,
    dBeta[iR,:],
    "o-",
    label="dMz/dbeta"
)

l2 = ax2.plot(
    tracks,
    dDelta[iR,:],
    "s-",
    label="dMz/ddelta"
)

ax1.set_xlabel("Front Track Width (in)")
ax1.set_ylabel("dMz/dbeta (ft-lb/rad)")
ax2.set_ylabel("dMz/ddelta (ft-lb/rad)")

ax1.set_title(
    f"Front Track Sensitivity - Rear Track = {tracks[iR]:.0f} in"
)

lines = l1+l2
ax1.legend(lines,[x.get_label() for x in lines])

ax1.grid(True)
plt.tight_layout()
plt.show()