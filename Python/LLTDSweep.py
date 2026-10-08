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

FrontDistribution_values = np.arange(0.30,0.71,0.01)

tire = Tire(vp.TireModel,vp.TirePressure_bar)

OriginalFrontRollStiffness = vp.FrontRollStiffness
OriginalRearRollStiffness = vp.RearRollStiffness

TotalRollStiffness = OriginalFrontRollStiffness+OriginalRearRollStiffness


# ===================================
# STORAGE
# ===================================

results = []
YMD_results = []


# ===================================
# PARAMETER SWEEP
# ===================================

for FrontDistribution in FrontDistribution_values:

    vp.FrontRollStiffness = TotalRollStiffness*FrontDistribution
    vp.RearRollStiffness = TotalRollStiffness*(1-FrontDistribution)

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

    dMz_dbeta = (Mz_grid[i0,j0+1]-Mz_grid[i0,j0-1])/(2*DEG2RAD)
    dMz_ddelta = (Mz_grid[i0+1,j0]-Mz_grid[i0-1,j0])/(2*DEG2RAD)

    CSR = dMz_ddelta/dMz_dbeta


    # ===================================
    # YMD LIMITS
    # ===================================

    Mz_ftlb_grid = Mz_grid*vp.NM2FTLB

    AyMax = np.max(Ay_grid)
    MzMax = np.max(Mz_ftlb_grid)

    AyMax_index = np.unravel_index(np.argmax(Ay_grid),Ay_grid.shape)
    MzAtAyMax = Mz_ftlb_grid[AyMax_index]


    # ===================================
    # SAVE RESULTS
    # ===================================

    results.append([
        FrontDistribution*100,
        AyMax,
        MzMax,
        MzAtAyMax,
        dMz_dbeta*vp.NM2FTLB,
        dMz_ddelta*vp.NM2FTLB,
        CSR
    ])

    YMD_results.append({
        "distribution":FrontDistribution*100,
        "Ay":Ay_grid.copy(),
        "Mz":Mz_ftlb_grid.copy()
    })

    print(
        f"Front Distribution = {FrontDistribution*100:.0f}% | "
        f"Ay Max = {AyMax:.3f}g | "
        f"CSR = {CSR:.2f}"
    )


# ===================================
# RESTORE VEHICLE PARAMETERS
# ===================================

vp.FrontRollStiffness = OriginalFrontRollStiffness
vp.RearRollStiffness = OriginalRearRollStiffness

results = np.array(results)


# ===================================
# PRINT RESULTS
# ===================================

print()
print("================================================================================================================")
print("ROLL STIFFNESS DISTRIBUTION SWEEP")
print("================================================================================================================")
print("Front %    Ay Max     Mz Max     Mz@AyMax     dMz/dbeta     dMz/ddelta     CSR")
print("----------------------------------------------------------------------------------------------------------------")

for row in results:

    print(
        f"{row[0]:6.1f}    "
        f"{row[1]:6.3f}    "
        f"{row[2]:8.1f}    "
        f"{row[3]:10.1f}    "
        f"{row[4]:11.1f}    "
        f"{row[5]:12.1f}    "
        f"{row[6]:6.2f}"
    )


# ===================================
# AY MAX PLOT
# ===================================

plt.figure(figsize=(10,7))

plt.plot(results[:,0],results[:,1],"o-")

plt.xlabel("Front Roll Stiffness Distribution (%)")
plt.ylabel("Maximum Lateral Acceleration (g)")
plt.title("Roll Stiffness Distribution vs Maximum Lateral Acceleration")

plt.grid(True)
plt.tight_layout()
plt.show()


# ===================================
# CSR PLOT
# ===================================

plt.figure(figsize=(10,7))

plt.plot(results[:,0],results[:,6],"o-")

plt.axhline(0,linewidth=1)

plt.xlabel("Front Roll Stiffness Distribution (%)")
plt.ylabel("Control / Stability Ratio")
plt.title("Roll Stiffness Distribution vs Control / Stability Ratio")

plt.grid(True)
plt.tight_layout()
plt.show()


# ===================================
# MZ AT AY MAX PLOT
# ===================================

plt.figure(figsize=(10,7))

plt.plot(results[:,0],results[:,3],"o-")

plt.axhline(0,linewidth=1)

plt.xlabel("Front Roll Stiffness Distribution (%)")
plt.ylabel("Yaw Moment at Maximum Ay (ft-lb)")
plt.title("Roll Stiffness Distribution vs Yaw Moment at Maximum Ay")

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

ax1.axhline(0,linewidth=1)

ax1.set_xlabel("Front Roll Stiffness Distribution (%)")
ax1.set_ylabel("dMz/dbeta (ft-lb/rad)")
ax2.set_ylabel("dMz/ddelta (ft-lb/rad)")

ax1.set_title("Roll Stiffness Distribution vs Yaw Moment Derivatives")

lines = line1+line2
labels = [line.get_label() for line in lines]

ax1.legend(lines,labels)

ax1.grid(True)

plt.tight_layout()
plt.show()


# ===================================
# INDIVIDUAL YMD COMPARISON
# ===================================

fig,axes = plt.subplots(3,3,figsize=(16,14))

AyMin = min(np.min(YMD["Ay"]) for YMD in YMD_results)
AyMax = max(np.max(YMD["Ay"]) for YMD in YMD_results)

MzMin = min(np.min(YMD["Mz"]) for YMD in YMD_results)
MzMax = max(np.max(YMD["Mz"]) for YMD in YMD_results)

for ax,YMD in zip(axes.flat,YMD_results):

    Ay = YMD["Ay"]
    Mz = YMD["Mz"]
    distribution = YMD["distribution"]

    for i in range(len(delta_values)):
        ax.plot(Ay[i,:],Mz[i,:],linewidth=1)

    ax.axhline(0,linewidth=0.8)
    ax.axvline(0,linewidth=0.8)

    ax.set_xlim(AyMin,AyMax)
    ax.set_ylim(MzMin,MzMax)

    ax.set_title(f"{distribution:.0f}% Front")

    ax.set_xlabel("Ay (g)")
    ax.set_ylabel("Mz (ft-lb)")

    ax.grid(True)

plt.suptitle(
    "YMD - Front Roll Stiffness Distribution",
    fontsize=16
)

plt.tight_layout()
plt.show()