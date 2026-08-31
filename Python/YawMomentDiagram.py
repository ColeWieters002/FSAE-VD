import numpy as np
import matplotlib.pyplot as plt
import VehicleParameters as vp
from VehicleParameters import DEG2RAD
from YMDSim import Tire,solve

def main():
    #Inputs
    Vx = 11.75
    beta_values=np.arange(-8,9,0.5)
    delta_values=np.arange(-8,9,0.5)

    tire=Tire(vp.TireModel,vp.TirePressure_bar)

    Ay_grid=np.zeros((len(delta_values),len(beta_values)))
    Mz_grid=np.zeros((len(delta_values),len(beta_values)))

    #Run YMD Sweep
    for i,delta_deg in enumerate(delta_values):
        for j,beta_deg in enumerate(beta_values):
            beta=beta_deg*DEG2RAD
            delta=delta_deg*DEG2RAD

            Ay,Mz=solve(Vx,beta,delta,vp,tire)

            Ay_grid[i,j]=Ay/vp.Gravity
            Mz_grid[i,j]=Mz

    #Plot YMD
    plt.figure(figsize=(11,8))

    #Constant Steering Angle Lines
    for i,delta_deg in enumerate(delta_values):
        plt.plot(Ay_grid[i,:],Mz_grid[i,:],label=f"δ={delta_deg}°")

    #Constant Sideslip Angle Lines
    for j,beta_deg in enumerate(beta_values):
        plt.plot(Ay_grid[:,j],Mz_grid[:,j],"--")

        if beta_deg%2==0:
            x=Ay_grid[-1,j]
            y=Mz_grid[-1,j]
            plt.text(x,y,f"β={beta_deg}°",fontsize=8)

    #Find Mz=0 Trim Points
    trim_Ay=[]
    trim_delta=[]

    for j,beta_deg in enumerate(beta_values):
        for i in range(len(delta_values)-1):
            Mz1=Mz_grid[i,j]
            Mz2=Mz_grid[i+1,j]

            if Mz1==0:
                trim_Ay.append(Ay_grid[i,j])
                trim_delta.append(delta_values[i])

            elif Mz1*Mz2<0:
                frac=-Mz1/(Mz2-Mz1)

                Ay_trim=Ay_grid[i,j]+frac*(Ay_grid[i+1,j]-Ay_grid[i,j])
                delta_trim=delta_values[i]+frac*(delta_values[i+1]-delta_values[i])

                trim_Ay.append(Ay_trim)
                trim_delta.append(delta_trim)

    #Plot Trim Points
    #if len(trim_Ay)>0:
        #plt.scatter(trim_Ay,np.zeros(len(trim_Ay)),s=35,zorder=5)

    plt.axhline(0,linewidth=1)
    plt.axvline(0,linewidth=1)

    plt.xlabel("Lateral Acceleration (g)")
    plt.ylabel("Yaw Moment (N·m)")
    plt.title(f"Yaw Moment Diagram - Vx={Vx:.2f} m/s")
    plt.grid(True)

    plt.legend(title="Constant Steering Angle",bbox_to_anchor=(1.02,1),loc="upper left")
    plt.tight_layout()
    plt.show()

    #Center Derivatives
    i0=np.where(delta_values==0)[0][0]
    j0=np.where(beta_values==0)[0][0]

    dMz_dbeta=(Mz_grid[i0,j0+1]-Mz_grid[i0,j0-1])/(2*DEG2RAD)
    dMz_ddelta=(Mz_grid[i0+1,j0]-Mz_grid[i0-1,j0])/(2*DEG2RAD)

    print()
    print("===================================")
    print("YMD CENTER DERIVATIVES")
    print("===================================")
    print(f"dMz/dbeta = {dMz_dbeta:.2f} N*m/rad")
    print(f"dMz/ddelta = {dMz_ddelta:.2f} N*m/rad")

    #Trim Steering Plot
    if len(trim_Ay)>0:
        trim_Ay=np.array(trim_Ay)
        trim_delta=np.array(trim_delta)

        order=np.argsort(trim_Ay)
        trim_Ay=trim_Ay[order]
        trim_delta=trim_delta[order]

        plt.figure(figsize=(8,6))
        plt.plot(trim_Ay,trim_delta,"o-")
        plt.axhline(0,linewidth=1)
        plt.axvline(0,linewidth=1)
        plt.xlabel("Lateral Acceleration (g)")
        plt.ylabel("Trim Steering Angle δ (deg)")
        plt.title("Trim Steering Angle vs Lateral Acceleration")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

if __name__=="__main__":
    main()