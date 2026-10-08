import numpy as np
import TireFunctions as MF
import VehicleParameters as vp
from VehicleParameters import LBF2N,N2LBF,FTLB2NM,NM2FTLB,FT2M,M2FT,IN2M,M2IN,RAD2DEG,DEG2RAD
from scipy.optimize import least_squares as ls
import matplotlib.pyplot as plt


class Tire:
    def __init__(self,tir_path,pressure_bar=None):
        sec=self._parse_tir(tir_path)

        self.p_fy=np.array([sec["LATERAL_COEFFICIENTS"][n] for n in MF.FY_Params],float)
        self.q_mz=np.array([sec["ALIGNING_COEFFICIENTS"][n] for n in MF.MZ_Params],float)
        self.p_fx=np.array([sec["LONGITUDINAL_COEFFICIENTS"][n] for n in MF.FX_Params],float)

        self.fz0=sec["VERTICAL"]["FNOMIN"]/LBF2N
        self.r0=sec["DIMENSION"]["UNLOADED_RADIUS"]/FT2M

        self.pressure_bar=pressure_bar
        self.path=tir_path


    def FY(self,SA_deg,FZ_lbf,IA_deg):
        return MF.FY((SA_deg,FZ_lbf,IA_deg),self.p_fy,self.fz0)


    def MZ(self,SA_deg,FZ_lbf,IA_deg):
        return MF.MZ((SA_deg,FZ_lbf,IA_deg),self.q_mz,self.p_fy,self.fz0,self.r0)


    def FX(self,SR,FZ_lbf,IA_deg):
        return MF.FX((SR,FZ_lbf,IA_deg),self.p_fx,self.fz0)


    def GraphSAvsFY(self,graphFz=100,graphInclinationAngle=10):
        graphSlipAngle=np.linspace(-20,20,50)
        graphLatForce=self.FY(graphSlipAngle,graphFz,graphInclinationAngle)

        plt.title("Slip Angle vs. Lateral Force")
        plt.plot(graphSlipAngle,graphLatForce)
        plt.grid(True)
        plt.show()


    @staticmethod
    def _parse_tir(path):
        sections,current={},None

        with open(path,"r") as f:
            for raw in f:
                line=raw.split("!",1)[0].strip()

                if not line:
                    continue

                if line.startswith("[") and line.endswith("]"):
                    current=line[1:-1].strip().upper()
                    sections[current]={}
                    continue

                if current is None or "=" not in line:
                    continue

                k,v=line.split("=",1)

                k=k.strip()
                v=v.strip().strip("'").strip('"')

                try:
                    v=float(v)
                except ValueError:
                    pass

                sections[current][k]=v

        return sections


def solve(Vx,delta,vp,tire,x0=None,max_iter=200,tol_beta=1e-4,tol_r=1e-3,tol_Ay=1e-3,tol_phi=1e-5,cost=1e-6,debug=False,return_data=False):

    #Base Variables
    L=vp.Wheelbase_mm/1000.0
    a=L*(1.0-vp.WeightDist)
    b=L*vp.WeightDist

    m=vp.TotalMass_kg
    ms=vp.SprungMass_kg
    g=vp.Gravity

    tf=vp.FTrackwidth_mm/1000.0
    tr=vp.RTrackwidth_mm/1000.0

    FrontRC=vp.FrontRollCenter_mm/1000.0
    RearRC=vp.RearRollCenter_mm/1000.0
    h_cg=vp.CG_mm/1000.0

    RollAxisHeight=(FrontRC*b+RearRC*a)/L

    Kphi_Total=vp.FrontRollStiffness+vp.RearRollStiffness


    def calculate_state(beta,r,Ay,phi,debug=False):

        #Slip Angles
        Vy=Vx*np.tan(beta)

        alpha_FL=delta-np.arctan2(Vy+r*a,Vx-r*tf/2)
        alpha_FR=delta-np.arctan2(Vy+r*a,Vx+r*tf/2)

        alpha_RL=-np.arctan2(Vy-r*b,Vx-r*tr/2)
        alpha_RR=-np.arctan2(Vy-r*b,Vx+r*tr/2)


        #Base FZs
        FZ_FL=m*g*vp.WeightDist/2
        FZ_FR=m*g*vp.WeightDist/2

        FZ_RL=m*g*(1.0-vp.WeightDist)/2
        FZ_RR=m*g*(1.0-vp.WeightDist)/2


        #Downforce
        DF=.5*vp.AirDensity*Vx**2*vp.CL*vp.A

        DF_FL=-DF*vp.AeroBalance*.5
        DF_FR=-DF*vp.AeroBalance*.5

        DF_RL=-DF*(1.0-vp.AeroBalance)*.5
        DF_RR=-DF*(1.0-vp.AeroBalance)*.5


        # Static camber is defined at static ride height.  In this
        # quasi-steady bicycle model, aero is the only symmetric axle-heave
        # displacement relative to that state.  Body roll is handled below
        # through the measured roll-camber gain, so do not add it here.
        if getattr(vp,"IncludeAeroHeaveCamber",True):
            FrontAeroDownforce=-DF*vp.AeroBalance
            RearAeroDownforce=-DF*(1.0-vp.AeroBalance)
            FrontAeroHeave_mm=1000.0*FrontAeroDownforce/vp.FrontHeaveStiffness
            RearAeroHeave_mm=1000.0*RearAeroDownforce/vp.RearHeaveStiffness
        else:
            FrontAeroHeave_mm=0.0
            RearAeroHeave_mm=0.0


        #Lateral Load Transfer
        Ys=ms*Ay

        YF=Ys*b/L
        YR=Ys*a/L

        FrontGeoMoment=YF*FrontRC
        RearGeoMoment=YR*RearRC

        FrontElasticMoment=vp.FrontRollStiffness*phi
        RearElasticMoment=vp.RearRollStiffness*phi

        FrontRollMoment=FrontGeoMoment+FrontElasticMoment
        RearRollMoment=RearGeoMoment+RearElasticMoment

        FrontLoadTransfer=FrontRollMoment/tf
        RearLoadTransfer=RearRollMoment/tr

        DF_FL-=FrontLoadTransfer
        DF_FR+=FrontLoadTransfer

        DF_RL-=RearLoadTransfer
        DF_RR+=RearLoadTransfer


        #Apply Aero and Load Transfer
        FZ_FL_raw=FZ_FL+DF_FL
        FZ_FR_raw=FZ_FR+DF_FR

        FZ_RL_raw=FZ_RL+DF_RL
        FZ_RR_raw=FZ_RR+DF_RR


        #Wheel Lift
        FZ_FL,FZ_FR=FZ_FL_raw,FZ_FR_raw

        if FZ_FL<0:
            FZ_FR+=FZ_FL
            FZ_FL=0.0

        elif FZ_FR<0:
            FZ_FL+=FZ_FR
            FZ_FR=0.0


        FZ_RL,FZ_RR=FZ_RL_raw,FZ_RR_raw

        if FZ_RL<0:
            FZ_RR+=FZ_RL
            FZ_RL=0.0

        elif FZ_RR<0:
            FZ_RL+=FZ_RR
            FZ_RR=0.0


        FZ_FL=max(FZ_FL,0.0)
        FZ_FR=max(FZ_FR,0.0)
        FZ_RL=max(FZ_RL,0.0)
        FZ_RR=max(FZ_RR,0.0)


        # Camber: static + aero heave + roll + steering gain.  The parameter
        # helper returns vehicle-coordinate camber; the right-tire TIR input
        # is mirrored below because the supplied TIR is a left-tire model.
        if not hasattr(vp,"CamberAtCorner_deg"):
            raise AttributeError(
                "VehicleParameters.py must define CamberAtCorner_deg. "
                "Use the matching VehicleParameters.py supplied with this upgrade."
            )

        roll_deg=phi*RAD2DEG
        steer_deg=delta*RAD2DEG
        gamma_FL=vp.CamberAtCorner_deg("front","left",roll_deg,steer_deg,FrontAeroHeave_mm)
        gamma_FR=vp.CamberAtCorner_deg("front","right",roll_deg,steer_deg,FrontAeroHeave_mm)
        gamma_RL=vp.CamberAtCorner_deg("rear","left",roll_deg,0.0,RearAeroHeave_mm)
        gamma_RR=vp.CamberAtCorner_deg("rear","right",roll_deg,0.0,RearAeroHeave_mm)

        # Physical wheel-local camber for reporting.  The right-side gamma
        # values above are sign-mirrored only for the left-tire TIR convention;
        # they are not a "positive camber" result at the right wheels.
        Camber_FL_deg=gamma_FL
        Camber_FR_deg=-gamma_FR
        Camber_RL_deg=gamma_RL
        Camber_RR_deg=-gamma_RR


        # Equivalent wheel-center travel relative to static ride height.
        # Positive is bump/compression; negative is rebound/extension.
        # This is rigid-body wheel travel, not damper stroke: convert through
        # the appropriate motion ratio before comparing to damper travel.
        # Positive roll in this model loads/compresses the right side.
        FrontRollTravel_mm=1000.0*(tf/2.0)*np.tan(phi)
        RearRollTravel_mm=1000.0*(tr/2.0)*np.tan(phi)
        Travel_FL_mm=FrontAeroHeave_mm-FrontRollTravel_mm
        Travel_FR_mm=FrontAeroHeave_mm+FrontRollTravel_mm
        Travel_RL_mm=RearAeroHeave_mm-RearRollTravel_mm
        Travel_RR_mm=RearAeroHeave_mm+RearRollTravel_mm


        #Tire Lateral Forces
        FY_FL=tire.FY(alpha_FL*RAD2DEG,FZ_FL*N2LBF,gamma_FL)*LBF2N
        FY_FR=-tire.FY(-alpha_FR*RAD2DEG,FZ_FR*N2LBF,-gamma_FR)*LBF2N

        FY_RL=tire.FY(alpha_RL*RAD2DEG,FZ_RL*N2LBF,gamma_RL)*LBF2N
        FY_RR=-tire.FY(-alpha_RR*RAD2DEG,FZ_RR*N2LBF,-gamma_RR)*LBF2N


        #Body Coordinates
        Fx_FL_body=-FY_FL*np.sin(delta)
        Fy_FL_body=FY_FL*np.cos(delta)
        Fx_FR_body=-FY_FR*np.sin(delta)
        Fy_FR_body=FY_FR*np.cos(delta)
        Fx_RL_body=0.0
        Fy_RL_body=FY_RL
        Fx_RR_body=0.0
        Fy_RR_body=FY_RR


        #Tire Aligning Moments
        TireMZ_FL=tire.MZ(alpha_FL*RAD2DEG,FZ_FL*N2LBF,gamma_FL)*FTLB2NM
        TireMZ_FR=-tire.MZ(-alpha_FR*RAD2DEG,FZ_FR*N2LBF,-gamma_FR)*FTLB2NM

        TireMZ_RL=tire.MZ(alpha_RL*RAD2DEG,FZ_RL*N2LBF,gamma_RL)*FTLB2NM
        TireMZ_RR=-tire.MZ(-alpha_RR*RAD2DEG,FZ_RR*N2LBF,-gamma_RR)*FTLB2NM

        TireMZ_Total=TireMZ_FL+TireMZ_FR+TireMZ_RL+TireMZ_RR


        #Yaw Moment From Tire Forces About CG
        ForceMoment_FL=a*Fy_FL_body-(tf/2)*Fx_FL_body
        ForceMoment_FR=a*Fy_FR_body-(-tf/2)*Fx_FR_body

        ForceMoment_RL=-b*Fy_RL_body
        ForceMoment_RR=-b*Fy_RR_body

        ForceMoment_Total=ForceMoment_FL+ForceMoment_FR+ForceMoment_RL+ForceMoment_RR


        #Totals
        FY_Total=Fy_FL_body+Fy_FR_body+Fy_RL_body+Fy_RR_body
        Mz_Total=ForceMoment_Total+TireMZ_Total

        LoadTransferTotal=FrontLoadTransfer+RearLoadTransfer

        if abs(LoadTransferTotal)>1e-12:
            TLLTD_Front=FrontLoadTransfer/LoadTransferTotal
            TLLTD_Rear=RearLoadTransfer/LoadTransferTotal
        else:
            TLLTD_Front=np.nan
            TLLTD_Rear=np.nan


        data={
            "Vx_mps":Vx,
            "delta_deg":delta*RAD2DEG,
            "beta_deg":beta*RAD2DEG,
            "r_rad_s":r,
            "r_deg_s":r*RAD2DEG,
            "Ay":Ay,
            "Ay_g":Ay/g,
            "phi_rad":phi,
            "phi_deg":phi*RAD2DEG,
            "Vy_mps":Vy,

            "FrontLoadTransfer_N":FrontLoadTransfer,
            "RearLoadTransfer_N":RearLoadTransfer,
            "TLLTD_Front":TLLTD_Front,
            "TLLTD_Rear":TLLTD_Rear,

            "RollAxisHeight_m":RollAxisHeight,
            "FrontGeoMoment_Nm":FrontGeoMoment,
            "RearGeoMoment_Nm":RearGeoMoment,
            "FrontElasticMoment_Nm":FrontElasticMoment,
            "RearElasticMoment_Nm":RearElasticMoment,
            "FrontRollMoment_Nm":FrontRollMoment,
            "RearRollMoment_Nm":RearRollMoment,

            "FrontAeroHeave_mm":FrontAeroHeave_mm,
            "RearAeroHeave_mm":RearAeroHeave_mm,
            "FrontRollTravel_mm":FrontRollTravel_mm,
            "RearRollTravel_mm":RearRollTravel_mm,
            "Travel_FL_mm":Travel_FL_mm,
            "Travel_FR_mm":Travel_FR_mm,
            "Travel_RL_mm":Travel_RL_mm,
            "Travel_RR_mm":Travel_RR_mm,

            "gamma_FL_deg":gamma_FL,
            "gamma_FR_deg":gamma_FR,
            "gamma_RL_deg":gamma_RL,
            "gamma_RR_deg":gamma_RR,
            "Camber_FL_deg":Camber_FL_deg,
            "Camber_FR_deg":Camber_FR_deg,
            "Camber_RL_deg":Camber_RL_deg,
            "Camber_RR_deg":Camber_RR_deg,

            "FZ_FL_N":FZ_FL,
            "FZ_FR_N":FZ_FR,
            "FZ_RL_N":FZ_RL,
            "FZ_RR_N":FZ_RR,

            "FZ_FL_raw_N":FZ_FL_raw,
            "FZ_FR_raw_N":FZ_FR_raw,
            "FZ_RL_raw_N":FZ_RL_raw,
            "FZ_RR_raw_N":FZ_RR_raw,

            "FY_FL_N":FY_FL,
            "FY_FR_N":FY_FR,
            "FY_RL_N":FY_RL,
            "FY_RR_N":FY_RR,
            "FY_Total_N":FY_Total,

            "alpha_FL_deg":alpha_FL*RAD2DEG,
            "alpha_FR_deg":alpha_FR*RAD2DEG,
            "alpha_RL_deg":alpha_RL*RAD2DEG,
            "alpha_RR_deg":alpha_RR*RAD2DEG,

            "TireMZ_FL_Nm":TireMZ_FL,
            "TireMZ_FR_Nm":TireMZ_FR,
            "TireMZ_RL_Nm":TireMZ_RL,
            "TireMZ_RR_Nm":TireMZ_RR,
            "TireMZ_Total_Nm":TireMZ_Total,

            "ForceMoment_Total_Nm":ForceMoment_Total,
            "Mz_Total_Nm":Mz_Total,
            "Mz_Total_ftlb":Mz_Total*NM2FTLB,
        }


        if debug:
            with open("Bicycle_Debug.txt","a") as f:
                f.write("\n====================================\n")
                f.write(f"Vx = {Vx}\n")
                f.write(f"delta = {delta*RAD2DEG} deg\n")
                f.write(f"beta = {beta*RAD2DEG} deg\n")
                f.write(f"r = {r} rad/s\n")
                f.write(f"Ay = {Ay}\n")
                f.write(f"Ay = {Ay/g} g\n")
                f.write(f"phi = {phi*RAD2DEG} deg\n")

                f.write("\nSlip Angles (deg)\n")
                f.write(f"FL = {data['alpha_FL_deg']}\n")
                f.write(f"FR = {data['alpha_FR_deg']}\n")
                f.write(f"RL = {data['alpha_RL_deg']}\n")
                f.write(f"RR = {data['alpha_RR_deg']}\n")

                f.write("\nCamber (deg; physical wheel-local, negative=negative camber)\n")
                f.write(f"FL = {Camber_FL_deg}\n")
                f.write(f"FR = {Camber_FR_deg}\n")
                f.write(f"RL = {Camber_RL_deg}\n")
                f.write(f"RR = {Camber_RR_deg}\n")

                f.write("\nWheel Travel (mm; +compression, -extension)\n")
                f.write(f"FL = {Travel_FL_mm}\n")
                f.write(f"FR = {Travel_FR_mm}\n")
                f.write(f"RL = {Travel_RL_mm}\n")
                f.write(f"RR = {Travel_RR_mm}\n")

                f.write("\nFZ (N)\n")
                f.write(f"FL = {FZ_FL}\n")
                f.write(f"FR = {FZ_FR}\n")
                f.write(f"RL = {FZ_RL}\n")
                f.write(f"RR = {FZ_RR}\n")

                f.write("\nFY (N)\n")
                f.write(f"FL = {FY_FL}\n")
                f.write(f"FR = {FY_FR}\n")
                f.write(f"RL = {FY_RL}\n")
                f.write(f"RR = {FY_RR}\n")
                f.write(f"FY Total = {FY_Total}\n")

                f.write("\nYaw Moments (Nm)\n")
                f.write(f"Force Moment Total = {ForceMoment_Total}\n")
                f.write(f"Tire MZ Total = {TireMZ_Total}\n")
                f.write(f"TOTAL Mz = {Mz_Total}\n")


        return FY_Total,Mz_Total,data


    def residual(x,delta):
        beta,r,Ay,phi=x

        FY_Total,Mz_Total,_=calculate_state(beta,r,Ay,phi,delta)

        #Lateral Force Equilibrium
        R1=(Ay-FY_Total/m)/g

        #Yaw Moment Equilibrium
        R2=Mz_Total/(m*g*L)

        #Steady State Kinematic Closure
        R3=(Ay-Vx*r)/g

        #Roll Equilibrium
        RollResidual=Kphi_Total*phi-ms*Ay*(h_cg-RollAxisHeight)
        R4=RollResidual/(m*g*L)

        return [R1,R2,R3,R4]


    #Initial Guess
    if x0 is None:
        x0=np.array([0.0,0.0,0.0,0.0])
        Step=2*DEG2RAD

        if delta>=0:
            DeltaSteps=np.arange(0,delta,Step)
        else:
            DeltaSteps=np.arange(0,delta,-Step)

        DeltaSteps=np.append(DeltaSteps,delta)

        for DeltaCurrent in DeltaSteps:
            result=ls(lambda x:residual(x,DeltaCurrent),x0,)
            x0=result.x
    else:
        x0=np.asarray(x0,float)

        if x0.shape!=(4,):
            raise ValueError("x0 must be [beta,r,Ay,phi]")


    #Solve beta, yaw rate, lateral acceleration, and roll angle
    result=ls(
        lambda x:residual(x,delta),
        x0,
        max_nfev=max_iter,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12
    )


    beta,r,Ay,phi=result.x
    FinalResidual=np.array(residual(result.x,delta))


    if (
        not result.success
        or result.cost>cost
        or abs(FinalResidual[0])>tol_Ay
        or abs(FinalResidual[1])>tol_r
        or abs(FinalResidual[2])>tol_Ay
        or abs(FinalResidual[3])>tol_phi
    ):
        raise RuntimeError(
            f"Solver failed to converge: cost={result.cost:.2e}, residual={FinalResidual}"
        )


    FY_Total,Mz_Total,data=calculate_state(beta,r,Ay,phi,debug=debug)

    data["solver_cost"]=result.cost
    data["solver_residuals"]=FinalResidual


    if return_data:
        return beta,r,Ay,phi,data

    return beta,r,Ay, data
