import numpy as np
import TireFunctions as MF
import VehicleParameters as vp
from VehicleParameters import LBF2N,N2LBF,FTLB2NM,NM2FTLB,FT2M,M2FT,IN2M,M2IN,RAD2DEG,DEG2RAD
from scipy.optimize import least_squares as ls


class Tire:
    def __init__(self, tir_path, pressure_bar=None):
        sec = self._parse_tir(tir_path)
        self.p_fy = np.array([sec["LATERAL_COEFFICIENTS"][n] for n in MF.FY_Params], float)
        self.q_mz = np.array([sec["ALIGNING_COEFFICIENTS"][n] for n in MF.MZ_Params], float)
        self.p_fx = np.array([sec["LONGITUDINAL_COEFFICIENTS"][n] for n in MF.FX_Params], float)
        self.fz0 = sec["VERTICAL"]["FNOMIN"] / LBF2N
        self.r0 = sec["DIMENSION"]["UNLOADED_RADIUS"] / FT2M
        self.pressure_bar = pressure_bar
        self.path = tir_path

    def FY(self, SA_deg, FZ_lbf, IA_deg):
        return MF.FY((SA_deg, FZ_lbf, IA_deg), self.p_fy, self.fz0)

    def MZ(self, SA_deg, FZ_lbf, IA_deg):
        return MF.MZ((SA_deg, FZ_lbf, IA_deg), self.q_mz, self.p_fy, self.fz0, self.r0)

    def FX(self, SR, FZ_lbf, IA_deg):
        return MF.FX((SR, FZ_lbf, IA_deg), self.p_fx, self.fz0)

    @staticmethod
    def _parse_tir(path):
        """Minimal .tir parser -> {SECTION: {KEY: value}} (same logic as TirePlotter)."""
        sections, current = {}, None
        with open(path, "r") as f:
            for raw in f:
                line = raw.split("!", 1)[0].strip()
                if not line:
                    continue
                if line.startswith("[") and line.endswith("]"):
                    current = line[1:-1].strip().upper()
                    sections[current] = {}
                    continue
                if current is None or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                k = k.strip()
                v = v.strip().strip("'").strip('"')
                try:
                    v = float(v)
                except ValueError:
                    pass
                sections[current][k] = v
        return sections


def solve(Vx, beta, delta, vp, tire):
    #Base Variables
    L = vp.Wheelbase_mm / 1000.0
    a = L * (1.0 - vp.WeightDist)
    b = L * vp.WeightDist
    m = vp.Mass_kg + vp.Driver_kg
    g = vp.Gravity

    tf = vp.FTrackwidth_mm / 1000.0
    tr = vp.RTrackwidth_mm / 1000.0

    FrontRC = vp.FrontRollCenter_mm / 1000.0
    RearRC = vp.RearRollCenter_mm / 1000.0
    h_cg = vp.CG_mm / 1000.0

    #Roll Stiffness Distribution
    Kphi_Total = vp.FrontRollStiffness + vp.RearRollStiffness
    FrontRollDistribution = vp.FrontRollStiffness / Kphi_Total
    RearRollDistribution = vp.RearRollStiffness / Kphi_Total

    def calculate_state(r, Ay):
        #Slip Angles
        Vy = Vx * np.tan(beta)

        alpha_FL = np.arctan2(Vy + r*a, Vx - r*tf/2) - delta
        alpha_FR = np.arctan2(Vy + r*a, Vx + r*tf/2) - delta
        alpha_RL = np.arctan2(Vy - r*b, Vx - r*tr/2)
        alpha_RR = np.arctan2(Vy - r*b, Vx + r*tr/2)

        #Base FZs (N)
        FZ_FL = m*g*vp.WeightDist/2
        FZ_FR = m*g*vp.WeightDist/2
        FZ_RL = m*g*(1-vp.WeightDist)/2
        FZ_RR = m*g*(1-vp.WeightDist)/2

        #Downforce per tire (N)
        DF_FL = -.5*vp.AirDensity*Vx**2*vp.CL*vp.A*vp.AeroBalance*.5
        DF_FR = -.5*vp.AirDensity*Vx**2*vp.CL*vp.A*vp.AeroBalance*.5
        DF_RL = -.5*vp.AirDensity*Vx**2*vp.CL*vp.A*(1-vp.AeroBalance)*.5
        DF_RR = -.5*vp.AirDensity*Vx**2*vp.CL*vp.A*(1-vp.AeroBalance)*.5

        #Lateral Load Transfer
        Y = m * Ay

        YF = Y * b / L
        YR = Y * a / L

        RollAxisHeight = (FrontRC*b + RearRC*a) / L

        FrontGeoMoment = YF * FrontRC
        RearGeoMoment = YR * RearRC

        ElasticMoment = Y * (h_cg - RollAxisHeight)

        FrontElasticMoment = ElasticMoment * FrontRollDistribution
        RearElasticMoment = ElasticMoment * RearRollDistribution

        FrontRollMoment = FrontGeoMoment + FrontElasticMoment
        RearRollMoment = RearGeoMoment + RearElasticMoment

        FrontLoadTransfer = FrontRollMoment / tf
        RearLoadTransfer = RearRollMoment / tr

        DF_FL -= FrontLoadTransfer
        DF_FR += FrontLoadTransfer
        DF_RL -= RearLoadTransfer
        DF_RR += RearLoadTransfer

        #Apply Aero and Load Transfer
        FZ_FL = max(FZ_FL + DF_FL, 0.0)
        FZ_FR = max(FZ_FR + DF_FR, 0.0)
        FZ_RL = max(FZ_RL + DF_RL, 0.0)
        FZ_RR = max(FZ_RR + DF_RR, 0.0)

        #Camber
        gamma_FL = vp.Camber_By_Travel_deg(0, "left")
        gamma_FR = vp.Camber_By_Travel_deg(0, "right")
        gamma_RL = vp.Camber_By_Travel_deg(0, "left")
        gamma_RR = vp.Camber_By_Travel_deg(0, "right")

        #Find FY
        FY_FL = tire.FY(
            alpha_FL*RAD2DEG,
            FZ_FL*N2LBF,
            gamma_FL
        )*LBF2N

        FY_FR = -tire.FY(
            -alpha_FR*RAD2DEG,
            FZ_FR*N2LBF,
            -gamma_FR
        )*LBF2N

        FY_RL = tire.FY(
            alpha_RL*RAD2DEG,
            FZ_RL*N2LBF,
            gamma_RL
        )*LBF2N

        FY_RR = -tire.FY(
            -alpha_RR*RAD2DEG,
            FZ_RR*N2LBF,
            -gamma_RR
        )*LBF2N

        #Find Tire Aligning Moments
        MZ_FL = tire.MZ(alpha_FL*RAD2DEG,FZ_FL*N2LBF,gamma_FL)*FTLB2NM

        MZ_FR = -tire.MZ(
            -alpha_FR*RAD2DEG,
            FZ_FR*N2LBF,
            -gamma_FR
        )*FTLB2NM

        MZ_RL = tire.MZ(alpha_RL*RAD2DEG,FZ_RL*N2LBF,gamma_RL)*FTLB2NM

        MZ_RR = -tire.MZ(
            -alpha_RR*RAD2DEG,
            FZ_RR*N2LBF,
            -gamma_RR
        )*FTLB2NM

        #Force and Moment Totals
        FYF = FY_FL + FY_FR
        FYR = FY_RL + FY_RR

        MZ_Tires = MZ_FL + MZ_FR + MZ_RL + MZ_RR

        FY_Total = FYF*np.cos(delta) + FYR

        Mz_Vehicle = a*FYF*np.cos(delta) - b*FYR + MZ_Tires

        print("alpha:")
        print(alpha_FL*RAD2DEG, alpha_FR*RAD2DEG, alpha_RL*RAD2DEG, alpha_RR*RAD2DEG)

        print("camber:")
        print(gamma_FL, gamma_FR, gamma_RL, gamma_RR)

        print("Fz:")
        print(FZ_FL, FZ_FR, FZ_RL, FZ_RR)

        print("Fy:")
        print(FY_FL, FY_FR, FY_RL, FY_RR)

        print("Mz tire:")
        print(MZ_FL, MZ_FR, MZ_RL, MZ_RR)

        return FY_Total, Mz_Vehicle

    def residual(x):
        r, Ay = x

        FY_Total, Mz_Vehicle = calculate_state(r, Ay)

        R1 = m*Ay - FY_Total
        R3 = Ay - Vx*r

        return [R1, R3]

    #Solve for Yaw Rate and Lateral Acceleration
    result = ls(residual, np.array([0.0, 0.0]))

    if result.cost > 1e-6:
        print(f"warning: residual not driven to zero (cost={result.cost:.2e})")

    r, Ay = result.x

    #Get Yaw Moment at Final Solved State
    FY_Total, Mz_Vehicle = calculate_state(r, Ay)

    

    return r, Ay, Mz_Vehicle