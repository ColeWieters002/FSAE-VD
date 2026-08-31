import numpy as np
from VehicleParameters import LBF2N,N2LBF,FTLB2NM,RAD2DEG

def calculate_state(Vx,delta,beta,r,Ay,vp,tire):
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

    #Tire Lateral Forces
    FY_FL=tire.FY(alpha_FL*RAD2DEG,FZ_FL*N2LBF,gamma_FL)*LBF2N
    FY_FR=-tire.FY(-alpha_FR*RAD2DEG,FZ_FR*N2LBF,-gamma_FR)*LBF2N
    FY_RL=tire.FY(alpha_RL*RAD2DEG,FZ_RL*N2LBF,gamma_RL)*LBF2N
    FY_RR=-tire.FY(-alpha_RR*RAD2DEG,FZ_RR*N2LBF,-gamma_RR)*LBF2N

    #Convert Front Tire Forces to Body Coordinates
    Fx_FL_body=-FY_FL*np.sin(delta)
    Fy_FL_body=FY_FL*np.cos(delta)
    Fx_FR_body=-FY_FR*np.sin(delta)
    Fy_FR_body=FY_FR*np.cos(delta)
    Fx_RL_body=0.0
    Fy_RL_body=FY_RL
    Fx_RR_body=0.0
    Fy_RR_body=FY_RR

    #Tire Self-Aligning Moments
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

    return {
        "alpha_FL":alpha_FL,
        "alpha_FR":alpha_FR,
        "alpha_RL":alpha_RL,
        "alpha_RR":alpha_RR,
        "gamma_FL":gamma_FL,
        "gamma_FR":gamma_FR,
        "gamma_RL":gamma_RL,
        "gamma_RR":gamma_RR,
        "FZ_FL":FZ_FL,
        "FZ_FR":FZ_FR,
        "FZ_RL":FZ_RL,
        "FZ_RR":FZ_RR,
        "FY_FL":FY_FL,
        "FY_FR":FY_FR,
        "FY_RL":FY_RL,
        "FY_RR":FY_RR,
        "Fx_FL_body":Fx_FL_body,
        "Fx_FR_body":Fx_FR_body,
        "Fx_RL_body":Fx_RL_body,
        "Fx_RR_body":Fx_RR_body,
        "Fy_FL_body":Fy_FL_body,
        "Fy_FR_body":Fy_FR_body,
        "Fy_RL_body":Fy_RL_body,
        "Fy_RR_body":Fy_RR_body,
        "TireMZ_FL":TireMZ_FL,
        "TireMZ_FR":TireMZ_FR,
        "TireMZ_RL":TireMZ_RL,
        "TireMZ_RR":TireMZ_RR,
        "ForceMoment_FL":ForceMoment_FL,
        "ForceMoment_FR":ForceMoment_FR,
        "ForceMoment_RL":ForceMoment_RL,
        "ForceMoment_RR":ForceMoment_RR,
        "FY_Total":FY_Total,
        "ForceMoment_Total":ForceMoment_Total,
        "TireMZ_Total":TireMZ_Total,
        "Mz_Total":Mz_Total
    }

def print_state(state):
    print("\n===================================")
    print("VEHICLE STATE")
    print("===================================")
    print(f"FL: alpha={state['alpha_FL']*RAD2DEG:.3f} deg, Fz={state['FZ_FL']:.1f} N, Fy={state['FY_FL']:.1f} N, Fy/Fz={abs(state['FY_FL'])/state['FZ_FL']:.3f}")
    print(f"FR: alpha={state['alpha_FR']*RAD2DEG:.3f} deg, Fz={state['FZ_FR']:.1f} N, Fy={state['FY_FR']:.1f} N, Fy/Fz={abs(state['FY_FR'])/state['FZ_FR']:.3f}")
    print(f"RL: alpha={state['alpha_RL']*RAD2DEG:.3f} deg, Fz={state['FZ_RL']:.1f} N, Fy={state['FY_RL']:.1f} N, Fy/Fz={abs(state['FY_RL'])/state['FZ_RL']:.3f}")
    print(f"RR: alpha={state['alpha_RR']*RAD2DEG:.3f} deg, Fz={state['FZ_RR']:.1f} N, Fy={state['FY_RR']:.1f} N, Fy/Fz={abs(state['FY_RR'])/state['FZ_RR']:.3f}")
    print(f"FY Total = {state['FY_Total']:.1f} N")
    print(f"Yaw Moment = {state['Mz_Total']:.3f} N*m")
