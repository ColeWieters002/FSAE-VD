import numpy as np
import matplotlib.pyplot as plt
import VehicleParameters as vp
from VehicleParameters import DEG2RAD
from YMDSim import Tire, solve


def main():

    # ===================================
    # INPUTS
    # ===================================

    Vx = 11.75

    beta_values = np.arange(-12,13,1)
    delta_values = np.arange(-12,13,1)

    open("YMD_Debug.txt","w").close()

    tire = Tire(vp.TireModel,vp.TirePressure_bar)


    # ===================================
    # CREATE STORAGE GRIDS
    # ===================================

    Ay_grid = np.zeros((len(delta_values),len(beta_values)))
    Mz_grid = np.zeros((len(delta_values),len(beta_values)))
    phi_grid = np.zeros((len(delta_values),len(beta_values)))


    # ===================================
    # RUN YMD SWEEP
    # ===================================

    for i,delta_deg in enumerate(delta_values):
        for j,beta_deg in enumerate(beta_values):

            beta = beta_deg*DEG2RAD
            delta = delta_deg*DEG2RAD

            Ay,Mz,phi,data = solve(Vx,beta,delta,vp,tire,True)

            Ay_grid[i,j] = Ay/vp.Gravity
            Mz_grid[i,j] = Mz
            phi_grid[i,j] = phi*vp.RAD2DEG


    # ===================================
    # CONVERT YAW MOMENT TO FT-LB
    # ===================================

    Mz_ftlb_grid = Mz_grid*vp.NM2FTLB


    # ===================================
    # LOCAL CONTROL / STABILITY DERIVATIVES
    # ===================================

    d_beta = (beta_values[1]-beta_values[0])*DEG2RAD
    d_delta = (delta_values[1]-delta_values[0])*DEG2RAD

    # axis=1 -> beta
    # axis=0 -> delta
    dMz_dbeta_grid = np.gradient(Mz_grid,d_beta,axis=1)
    dMz_ddelta_grid = np.gradient(Mz_grid,d_delta,axis=0)


    # ===================================
    # FIND CENTER POINT
    # ===================================

    i0 = np.where(delta_values==0)[0][0]
    j0 = np.where(beta_values==0)[0][0]


    # ===================================
    # CONTROL / STABILITY RATIO
    # ===================================

    with np.errstate(divide="ignore",invalid="ignore"):
        CSR_grid = dMz_ddelta_grid/dMz_dbeta_grid

    CSR_limit = 15.0

    CSR_plot = np.nan_to_num(
        CSR_grid,
        nan=0.0,
        posinf=CSR_limit,
        neginf=-CSR_limit
    )

    CSR_plot = np.clip(CSR_plot,-CSR_limit,CSR_limit)


    # ===================================
    # 3D YMD - CONTROL / STABILITY RATIO
    # ===================================

    fig = plt.figure(figsize=(14,10))
    ax = fig.add_subplot(111,projection="3d")

    X = Ay_grid
    Y = Mz_ftlb_grid
    Z = CSR_plot

    surf = ax.plot_surface(
        X,Y,Z,
        cmap="viridis",
        vmin=-CSR_limit,
        vmax=CSR_limit,
        edgecolor="black",
        linewidth=0.25,
        alpha=0.90,
        rstride=1,
        cstride=1
    )

    ax.set_xlabel("Lateral Acceleration $A_y$ (g)",labelpad=12)
    ax.set_ylabel("Yaw Moment $M_z$ (ft-lb)",labelpad=12)
    ax.set_zlabel("Control / Stability Ratio",labelpad=12)
    ax.set_title(f"3D Yaw Moment Diagram - Vx={Vx:.2f} m/s")

    fig.colorbar(
        surf,
        ax=ax,
        shrink=0.65,
        pad=0.10,
        label="Control / Stability Ratio"
    )

    ax.view_init(elev=90,azim=-90)
    ax.grid(True)

    plt.tight_layout()
    plt.show()


    # ===================================
    # CENTER DERIVATIVES
    # ===================================

    dMz_dbeta = (Mz_grid[i0,j0+1]-Mz_grid[i0,j0-1])/(2*DEG2RAD)
    dMz_ddelta = (Mz_grid[i0+1,j0]-Mz_grid[i0-1,j0])/(2*DEG2RAD)

    CSR_center = dMz_ddelta/dMz_dbeta


    # ===================================
    # PRINT CENTER DERIVATIVES
    # ===================================

    print()
    print("===================================")
    print("YMD CENTER DERIVATIVES")
    print("===================================")
    print(f"dMz/dbeta = {dMz_dbeta*vp.NM2FTLB:.2f} ft-lbs/rad")
    print(f"dMz/ddelta = {dMz_ddelta*vp.NM2FTLB:.2f} ft-lbs/rad")
    print(f"Control / Stability Ratio = {CSR_center:.2f}")


    # ===================================
    # PRINT YMD LIMITS
    # ===================================

    print()
    print("===================================")
    print("YMD MZ AND AY LIMITS")
    print("===================================")
    print(f"MZ MAX = {np.max(Mz_ftlb_grid):.2f} ft-lbs")
    print(f"MZ MIN = {np.min(Mz_ftlb_grid):.2f} ft-lbs")
    print(f"AY MAX = {np.max(Ay_grid):.2f} g's")
    print(f"AY MIN = {np.min(Ay_grid):.2f} g's")


    # ===================================
    # FIND Mz = 0 TRIM POINTS
    # ===================================

    trim_Ay = []
    trim_delta = []

    for j,beta_deg in enumerate(beta_values):
        for i in range(len(delta_values)-1):

            Mz1 = Mz_grid[i,j]
            Mz2 = Mz_grid[i+1,j]

            # Exact zero
            if Mz1==0:
                trim_Ay.append(Ay_grid[i,j])
                trim_delta.append(delta_values[i])

            # Zero crossing
            elif Mz1*Mz2<0:

                frac = -Mz1/(Mz2-Mz1)

                Ay_trim = Ay_grid[i,j]+frac*(Ay_grid[i+1,j]-Ay_grid[i,j])
                delta_trim = delta_values[i]+frac*(delta_values[i+1]-delta_values[i])

                trim_Ay.append(Ay_trim)
                trim_delta.append(delta_trim)


    # ===================================
    # TRIM STEERING PLOT
    # ===================================

    if len(trim_Ay)>0:

        trim_Ay = np.array(trim_Ay)
        trim_delta = np.array(trim_delta)

        order = np.argsort(trim_Ay)

        trim_Ay = trim_Ay[order]
        trim_delta = trim_delta[order]

        plt.figure(figsize=(10,7))

        plt.plot(trim_Ay,trim_delta,"o-")

        plt.axhline(0,linewidth=1)
        plt.axvline(0,linewidth=1)

        plt.xlabel("Lateral Acceleration (g)")
        plt.ylabel("Trim Steering Angle δ (deg)")
        plt.title("Trim Steering Angle vs Lateral Acceleration")

        plt.grid(True)
        plt.tight_layout()
        plt.show()


    # ===================================
    # SOLID STABILITY DERIVATIVE MAP
    # ===================================

    fig,ax = plt.subplots(figsize=(11,8))

    stability_map = ax.pcolormesh(
        Ay_grid,
        Mz_ftlb_grid,
        dMz_dbeta_grid,
        shading="auto",
        cmap="viridis"
    )

    fig.colorbar(
        stability_map,
        ax=ax,
        label=r"$dM_z/d\beta$ (N·m/rad)"
    )

    ax.axhline(0,linewidth=1)
    ax.axvline(0,linewidth=1)

    ax.set_xlabel("Lateral Acceleration $A_y$ (g)")
    ax.set_ylabel("Yaw Moment $M_z$ (ft-lb)")
    ax.set_title("Local Stability Derivative")

    ax.grid(True)

    plt.tight_layout()
    plt.show()


    # ===================================
    # SOLID CONTROL DERIVATIVE MAP
    # ===================================

    fig,ax = plt.subplots(figsize=(11,8))

    control_map = ax.pcolormesh(
        Ay_grid,
        Mz_ftlb_grid,
        dMz_ddelta_grid,
        shading="auto",
        cmap="viridis"
    )

    fig.colorbar(
        control_map,
        ax=ax,
        label=r"$dM_z/d\delta$ (N·m/rad)"
    )

    ax.axhline(0,linewidth=1)
    ax.axvline(0,linewidth=1)

    ax.set_xlabel("Lateral Acceleration $A_y$ (g)")
    ax.set_ylabel("Yaw Moment $M_z$ (ft-lb)")
    ax.set_title("Local Control Derivative")

    ax.grid(True)

    plt.tight_layout()
    plt.show()


if __name__=="__main__":
    main()