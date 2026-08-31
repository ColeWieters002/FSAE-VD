import numpy as np
import TireFunctions as MF
import VehicleParameters as vp
from VehicleParameters import LBF2N,N2LBF,FTLB2NM,NM2FTLB,FT2M,M2FT,IN2M,M2IN,RAD2DEG,DEG2RAD
from scipy.optimize import least_squares as ls

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

def solve(Vx,beta,delta,vp,tire,debug=False):
    #Base Variables
    L=vp.Wheelbase_mm/1000.0
    a=L*(1.0-vp.WeightDist)
    b=L*vp.WeightDist
    m=vp.Mass_kg+vp.Driver_kg
    g=vp.Gravity
    tf=vp.FTrackwidth_mm/1000.0
    tr=vp.RTrackwidth_mm/1000.0
    FrontRC=vp.FrontRollCenter_mm/1000.0
    RearRC=vp.RearRollCenter_mm/1000.0
    h_cg=vp.CG_mm/1000.0

    #Roll Stiffness Distribution
    Kphi_Total=vp.FrontRollStiffness+vp.RearRollStiffness
    FrontRollDistribution=vp.FrontRollStiffness/Kphi_Total
    RearRollDistribution=vp.RearRollStiffness/Kphi_Total

    def calculate_state(Ay,debug=False):
        #Zero Yaw Rate YMD
        r=0.0

        #Slip Angles
        Vy=Vx*np.tan(beta)

        alpha_FL=delta-np.arctan2(Vy+r*a,Vx-r*tf/2)
        alpha_FR=delta-np.arctan2(Vy+r*a,Vx+r*tf/2)
        alpha_RL=-np.arctan2(Vy-r*b,Vx-r*tr/2)
        alpha_RR=-np.arctan2(Vy-r*b,Vx+r*tr/2)

        #Base FZs
        FZ_FL=m*g*vp.WeightDist/2
        FZ_FR=m*g*vp.WeightDist/2
        FZ_RL=m*g*(1-vp.WeightDist)/2
        FZ_RR=m*g*(1-vp.WeightDist)/2

        #Downforce
        DF_FL=-.5*vp.AirDensity*Vx**2*vp.CL*vp.A*vp.AeroBalance*.5
        DF_FR=-.5*vp.AirDensity*Vx**2*vp.CL*vp.A*vp.AeroBalance*.5
        DF_RL=-.5*vp.AirDensity*Vx**2*vp.CL*vp.A*(1-vp.AeroBalance)*.5
        DF_RR=-.5*vp.AirDensity*Vx**2*vp.CL*vp.A*(1-vp.AeroBalance)*.5

        #Lateral Load Transfer
        Y=m*Ay
        YF=Y*b/L
        YR=Y*a/L

        RollAxisHeight=(FrontRC*b+RearRC*a)/L
        FrontGeoMoment=YF*FrontRC
        RearGeoMoment=YR*RearRC
        ElasticMoment=Y*(h_cg-RollAxisHeight)
        FrontElasticMoment=ElasticMoment*FrontRollDistribution
        RearElasticMoment=ElasticMoment*RearRollDistribution
        FrontRollMoment=FrontGeoMoment+FrontElasticMoment
        RearRollMoment=RearGeoMoment+RearElasticMoment
        FrontLoadTransfer=FrontRollMoment/tf
        RearLoadTransfer=RearRollMoment/tr

        DF_FL-=FrontLoadTransfer
        DF_FR+=FrontLoadTransfer
        DF_RL-=RearLoadTransfer
        DF_RR+=RearLoadTransfer

        #Apply Aero and Load Transfer
        FZ_FL=max(FZ_FL+DF_FL,0.0)
        FZ_FR=max(FZ_FR+DF_FR,0.0)
        FZ_RL=max(FZ_RL+DF_RL,0.0)
        FZ_RR=max(FZ_RR+DF_RR,0.0)

        #Camber
        gamma_FL=vp.Camber_By_Travel_deg(0,"left")
        gamma_FR=vp.Camber_By_Travel_deg(0,"right")
        gamma_RL=vp.Camber_By_Travel_deg(0,"left")
        gamma_RR=vp.Camber_By_Travel_deg(0,"right")

        #Find Tire Lateral Forces
        FY_FL=tire.FY(alpha_FL*RAD2DEG,FZ_FL*N2LBF,gamma_FL)*LBF2N
        FY_FR=-tire.FY(-alpha_FR*RAD2DEG,FZ_FR*N2LBF,-gamma_FR)*LBF2N
        FY_RL=tire.FY(alpha_RL*RAD2DEG,FZ_RL*N2LBF,gamma_RL)*LBF2N
        FY_RR=-tire.FY(-alpha_RR*RAD2DEG,FZ_RR*N2LBF,-gamma_RR)*LBF2N

        #Convert Tire Forces to Body Coordinates
        Fx_FL_body=-FY_FL*np.sin(delta)
        Fy_FL_body=FY_FL*np.cos(delta)
        Fx_FR_body=-FY_FR*np.sin(delta)
        Fy_FR_body=FY_FR*np.cos(delta)
        Fx_RL_body=0.0
        Fy_RL_body=FY_RL
        Fx_RR_body=0.0
        Fy_RR_body=FY_RR

        #Find Tire Self-Aligning Moments
        TireMZ_FL=tire.MZ(alpha_FL*RAD2DEG,FZ_FL*N2LBF,gamma_FL)*FTLB2NM
        TireMZ_FR=-tire.MZ(-alpha_FR*RAD2DEG,FZ_FR*N2LBF,-gamma_FR)*FTLB2NM
        TireMZ_RL=tire.MZ(alpha_RL*RAD2DEG,FZ_RL*N2LBF,gamma_RL)*FTLB2NM
        TireMZ_RR=-tire.MZ(-alpha_RR*RAD2DEG,FZ_RR*N2LBF,-gamma_RR)*FTLB2NM

        TireMZ_Total=TireMZ_FL+TireMZ_FR+TireMZ_RL+TireMZ_RR

        #Find Yaw Moment From Tire Forces About CG
        ForceMoment_FL=a*Fy_FL_body-(tf/2)*Fx_FL_body
        ForceMoment_FR=a*Fy_FR_body-(-tf/2)*Fx_FR_body
        ForceMoment_RL=-b*Fy_RL_body
        ForceMoment_RR=-b*Fy_RR_body

        ForceMoment_Total=ForceMoment_FL+ForceMoment_FR+ForceMoment_RL+ForceMoment_RR

        #Force and Moment Totals
        FY_Total=Fy_FL_body+Fy_FR_body+Fy_RL_body+Fy_RR_body
        Mz_Vehicle=ForceMoment_Total+TireMZ_Total

        if debug:
            with open("YMD_Debug.txt","a") as f:
                f.write("\n====================================\n")
                f.write(f"Vx = {Vx}\n")
                f.write(f"beta = {beta*RAD2DEG} deg\n")
                f.write(f"delta = {delta*RAD2DEG} deg\n")
                f.write("--- STATE ---\n")
                f.write(f"r = {r}\n")
                f.write(f"Ay = {Ay}\n")

                f.write("\nSlip Angles (deg)\n")
                f.write(f"FL = {alpha_FL*RAD2DEG}\n")
                f.write(f"FR = {alpha_FR*RAD2DEG}\n")
                f.write(f"RL = {alpha_RL*RAD2DEG}\n")
                f.write(f"RR = {alpha_RR*RAD2DEG}\n")

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

                f.write("\nForce Moments About CG (Nm)\n")
                f.write(f"FL = {ForceMoment_FL}\n")
                f.write(f"FR = {ForceMoment_FR}\n")
                f.write(f"RL = {ForceMoment_RL}\n")
                f.write(f"RR = {ForceMoment_RR}\n")
                f.write(f"Force Moment Total = {ForceMoment_Total}\n")

                f.write("\nTire Aligning Moments (Nm)\n")
                f.write(f"FL = {TireMZ_FL}\n")
                f.write(f"FR = {TireMZ_FR}\n")
                f.write(f"RL = {TireMZ_RL}\n")
                f.write(f"RR = {TireMZ_RR}\n")
                f.write(f"Tire MZ Total = {TireMZ_Total}\n")

                f.write(f"\nTOTAL Mz = {Mz_Vehicle}\n")

        return FY_Total,Mz_Vehicle

    def residual(x):
        Ay=x[0]
        FY_Total,Mz_Vehicle=calculate_state(Ay)
        R1=Ay-FY_Total/m
        return [R1]

    #Solve for Lateral Acceleration
    result=ls(residual,np.array([0.0]))

    if result.cost>1e-6:
        print(f"warning: residual not driven to zero (cost={result.cost:.2e})")

    Ay=result.x[0]

    #Get Yaw Moment at Final Solved State
    FY_Total,Mz_Vehicle=calculate_state(Ay,debug=debug)

    return Ay,Mz_Vehicle