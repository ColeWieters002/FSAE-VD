import numpy as np
import matplotlib.pyplot as plt
import VehicleParameters as vp
from VehicleParameters import DEG2RAD
from YMDSim import Tire,solve


# ===================================
# INPUTS
# ===================================

Vx = 11.75

beta_values = np.arange(-12,13,1)
delta_values = np.arange(-12,13,1)

CG_values_in = np.arange(7,15.01,0.25)
YMD_CG_values = [7,8,9,10,11,12,13,14,15]

FineStep = 0.1

tire = Tire(vp.TireModel,vp.TirePressure_bar)

OriginalCG = vp.CG_mm


# ===================================
# STORAGE
# ===================================

results = []
YMD_results = []


# ===================================
# PARAMETER SWEEP
# ===================================

for CG_in in CG_values_in:

    vp.CG_mm = CG_in*vp.IN2M*1000

    Ay_grid = np.zeros((len(delta_values),len(beta_values)))
    Mz_grid = np.zeros((len(delta_values),len(beta_values)))

    for i,delta_deg in enumerate(delta_values):
        for j,beta_deg in enumerate(beta_values):

            beta = beta_deg*DEG2RAD
            delta = delta_deg*DEG2RAD

            Ay,Mz,phi,data = solve(Vx,beta,delta,vp,tire,False)

            Ay_grid[i,j] = Ay/vp.Gravity
            Mz_grid[i,j] = Mz


    # ===================================
    # CENTER DERIVATIVES
    # ===================================

    i0 = np.where(delta_values==0)[0][0]
    j0 = np.where(beta_values==0)[0][0]

    dMz_dbeta = (
        Mz_grid[i0,j0+1]-Mz_grid[i0,j0-1]
    )/(2*DEG2RAD)

    dMz_ddelta = (
        Mz_grid[i0+1,j0]-Mz_grid[i0-1,j0]
    )/(2*DEG2RAD)

    CSR = dMz_ddelta/dMz_dbeta


    # ===================================
    # YMD LIMITS
    # ===================================

    Mz_ftlb_grid = Mz_grid*vp.NM2FTLB

    MzMax = np.max(Mz_ftlb_grid)


    # ===================================
    # COARSE AY MAX
    # ===================================

    i_max,j_max = np.unravel_index(
        np.argmax(Ay_grid),
        Ay_grid.shape
    )

    CoarseBeta = beta_values[j_max]
    CoarseDelta = delta_values[i_max]


    # ===================================
    # FINE AY MAX SEARCH
    # ===================================

    beta_fine = np.arange(
        max(-12,CoarseBeta-1),
        min(12,CoarseBeta+1)+FineStep/2,
        FineStep
    )

    delta_fine = np.arange(
        max(-12,CoarseDelta-1),
        min(12,CoarseDelta+1)+FineStep/2,
        FineStep
    )

    AyMax = -np.inf
    MzAtAyMax = 0.0
    BetaAtAyMax = 0.0
    DeltaAtAyMax = 0.0

    for delta_deg in delta_fine:
        for beta_deg in beta_fine:

            beta = beta_deg*DEG2RAD
            delta = delta_deg*DEG2RAD

            Ay,Mz,phi,data = solve(
                Vx,
                beta,
                delta,
                vp,
                tire,
                False
            )

            Ay_g = Ay/vp.Gravity

            if Ay_g > AyMax:

                AyMax = Ay_g
                MzAtAyMax = Mz*vp.NM2FTLB

                BetaAtAyMax = beta_deg
                DeltaAtAyMax = delta_deg


    # ===================================
    # SAVE RESULTS
    # ===================================

    results.append([
        CG_in,
        AyMax,
        MzMax,
        MzAtAyMax,
        dMz_dbeta*vp.NM2FTLB,
        dMz_ddelta*vp.NM2FTLB,
        CSR,
        BetaAtAyMax,
        DeltaAtAyMax
    ])

    if any(np.isclose(CG_in,value) for value in YMD_CG_values):

        YMD_results.append({
            "CG":CG_in,
            "Ay":Ay_grid.copy(),
            "Mz":Mz_ftlb_grid.copy()
        })

    print(
        f"CG = {CG_in:.2f} in | "
        f"Ay Max = {AyMax:.4f}g | "
        f"Mz = {MzAtAyMax:.1f} ft-lb | "
        f"Beta = {BetaAtAyMax:.1f} deg | "
        f"Delta = {DeltaAtAyMax:.1f} deg | "
        f"CSR = {CSR:.2f}"
    )


# ===================================
# RESTORE VEHICLE PARAMETERS
# ===================================

vp.CG_mm = OriginalCG

results = np.array(results)


# ===================================
# PRINT RESULTS
# ===================================

print()
print("================================================================================================================================")
print("CG HEIGHT SWEEP")
print("================================================================================================================================")
print("CG (in)   Ay Max    Mz Max   Mz@AyMax   dMz/dbeta   dMz/ddelta    CSR     Beta     Delta")
print("--------------------------------------------------------------------------------------------------------------------------------")

for row in results:

    print(
        f"{row[0]:7.2f}   "
        f"{row[1]:6.4f}   "
        f"{row[2]:7.1f}   "
        f"{row[3]:9.1f}   "
        f"{row[4]:10.1f}   "
        f"{row[5]:11.1f}   "
        f"{row[6]:6.2f}   "
        f"{row[7]:6.1f}   "
        f"{row[8]:6.1f}"
    )


# ===================================
# AY MAX PLOT
# ===================================

plt.figure(figsize=(10,7))

plt.plot(
    results[:,0],
    results[:,1],
    "o-"
)

plt.axvline(
    OriginalCG/1000/vp.IN2M,
    linewidth=1,
    linestyle="--",
    label="Current CG"
)

plt.xlabel("CG Height (in)")
plt.ylabel("Maximum Lateral Acceleration (g)")
plt.title("CG Height vs Maximum Lateral Acceleration")

plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ===================================
# CSR PLOT
# ===================================

plt.figure(figsize=(10,7))

plt.plot(
    results[:,0],
    results[:,6],
    "o-"
)

plt.axhline(
    0,
    linewidth=1
)

plt.axvline(
    OriginalCG/1000/vp.IN2M,
    linewidth=1,
    linestyle="--",
    label="Current CG"
)

plt.xlabel("CG Height (in)")
plt.ylabel("Control / Stability Ratio")
plt.title("CG Height vs Control / Stability Ratio")

plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ===================================
# MZ AT AY MAX PLOT
# ===================================

plt.figure(figsize=(10,7))

plt.plot(
    results[:,0],
    results[:,3],
    "o-"
)

plt.axhline(
    0,
    linewidth=1
)

plt.axvline(
    OriginalCG/1000/vp.IN2M,
    linewidth=1,
    linestyle="--",
    label="Current CG"
)

plt.xlabel("CG Height (in)")
plt.ylabel("Yaw Moment at Maximum Ay (ft-lb)")
plt.title("CG Height vs Yaw Moment at Maximum Ay")

plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ===================================
# CONTROL / STABILITY DERIVATIVES
# ===================================

fig,ax1 = plt.subplots(figsize=(10,7))
ax2 = ax1.twinx()

line1 = ax1.plot(
    results[:,0],
    results[:,4],
    "o-",
    label="dMz/dbeta"
)

line2 = ax2.plot(
    results[:,0],
    results[:,5],
    "s-",
    label="dMz/ddelta"
)

ax1.axvline(
    OriginalCG/1000/vp.IN2M,
    linewidth=1,
    linestyle="--",
    label="Current CG"
)

ax1.set_xlabel("CG Height (in)")
ax1.set_ylabel("dMz/dbeta (ft-lb/rad)")
ax2.set_ylabel("dMz/ddelta (ft-lb/rad)")

ax1.set_title("CG Height vs Yaw Moment Derivatives")


# ===================================
# MATCH RELATIVE Y-AXIS SCALES
# ===================================

beta_min = np.min(results[:,4])
beta_max = np.max(results[:,4])

delta_min = np.min(results[:,5])
delta_max = np.max(results[:,5])

beta_range = beta_max-beta_min
delta_range = delta_max-delta_min

beta_padding = beta_range*0.15
delta_padding = delta_range*0.15

ax1.set_ylim(
    beta_min-beta_padding,
    beta_max+beta_padding
)

ax2.set_ylim(
    delta_min-delta_padding,
    delta_max+delta_padding
)


# ===================================
# LEGEND
# ===================================

lines = line1+line2
labels = [line.get_label() for line in lines]

ax1.legend(
    lines,
    labels
)

ax1.grid(True)

plt.tight_layout()
plt.show()


# ===================================
# INDIVIDUAL YMD COMPARISON
# ===================================

fig,axes = plt.subplots(
    3,
    3,
    figsize=(16,14)
)

AyMin = min(
    np.min(YMD["Ay"])
    for YMD in YMD_results
)

AyMaxPlot = max(
    np.max(YMD["Ay"])
    for YMD in YMD_results
)

MzMin = min(
    np.min(YMD["Mz"])
    for YMD in YMD_results
)

MzMaxPlot = max(
    np.max(YMD["Mz"])
    for YMD in YMD_results
)

for ax,YMD in zip(axes.flat,YMD_results):

    Ay = YMD["Ay"]
    Mz = YMD["Mz"]
    CG = YMD["CG"]

    for i in range(len(delta_values)):

        ax.plot(
            Ay[i,:],
            Mz[i,:],
            linewidth=1
        )

    ax.axhline(
        0,
        linewidth=0.8
    )

    ax.axvline(
        0,
        linewidth=0.8
    )

    ax.set_xlim(
        AyMin,
        AyMaxPlot
    )

    ax.set_ylim(
        MzMin,
        MzMaxPlot
    )

    ax.set_title(
        f"CG Height = {CG:.0f} in"
    )

    ax.set_xlabel("Ay (g)")
    ax.set_ylabel("Mz (ft-lb)")

    ax.grid(True)

plt.suptitle(
    "YMD - CG Height Sweep",
    fontsize=16
)

plt.tight_layout()
plt.show()