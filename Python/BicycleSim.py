import numpy as np
import TireFunctions as MF
import VehicleParameters as vp
from scipy.optimize import least_squares as ls
import matplotlib.pyplot as plt

# Unit Conversions
LBF2N   = 4.4482216
N2LBF   = 1.0 / LBF2N
FTLB2NM = 1.3558
FT2M    = 0.3048
RAD2DEG = 180.0 / np.pi
DEG2RAD = np.pi / 180.0

class Tire:
    def __init__(self, tir_path, pressure_bar=None):
        sec = self._parse_tir(tir_path)
        self.p_fy = np.array([sec["LATERAL_COEFFICIENTS"][n]      for n in MF.FY_Params], float)
        self.q_mz = np.array([sec["ALIGNING_COEFFICIENTS"][n]     for n in MF.MZ_Params], float)
        self.p_fx = np.array([sec["LONGITUDINAL_COEFFICIENTS"][n] for n in MF.FX_Params], float)
        self.fz0  = sec["VERTICAL"]["FNOMIN"] / LBF2N          # nominal load [lbf]
        self.r0   = sec["DIMENSION"]["UNLOADED_RADIUS"] / FT2M  # radius [ft]
        self.pressure_bar = pressure_bar
        self.path = tir_path

    def FY(self, SA_deg, FZ_lbf, IA_deg):
        return MF.FY((SA_deg, FZ_lbf, IA_deg), self.p_fy, self.fz0)

    def MZ(self, SA_deg, FZ_lbf, IA_deg):
        return MF.MZ((SA_deg, FZ_lbf, IA_deg), self.q_mz, self.p_fy, self.fz0, self.r0)

    def FX(self, SR, FZ_lbf, IA_deg):
        return MF.FX((SR, FZ_lbf, IA_deg), self.p_fx, self.fz0)
    
    def GraphSAvsFY(self, graphFz = 100, graphInclinationAngle = 10):
        graphFz *= LBF2N
        graphInclinationAngle *= DEG2RAD
        graphSlipAngle = np.array(np.linspace(-20, 20, 50))
        graphLatForce = self.FY(graphSlipAngle,graphFz,graphInclinationAngle)
        plt.title("Slip Angle vs. Lateral Force")
        plt.plot(graphSlipAngle,graphLatForce)
        plt.show()

    @staticmethod
    def _parse_tir(path):
        """Minimal .tir parser -> {SECTION: {KEY: value}} (same logic as TirePlotter)."""
        sections, current = {}, None
        with open(path, "r") as f:
            for raw in f:
                line = raw.split("!", 1)[0].strip()      # drop comments
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


def solve(Vx, delta, vp, tire, max_iter=100, relax=0.4, tol_beta=1e-4, tol_r=1e-3, tol_Ay=1e-3):
    #Base Variables
    L = vp.Wheelbase_mm / 1000.0
    a = L * (1.0 - vp.WeightDist)     #CG to Front Axle (m)
    b = L * vp.WeightDist             #CG to Rear Axle (m)
    m = vp.Mass_kg + vp.Driver_kg     #Total mass (kg)
    g = vp.Gravity

    def residual(x):
        #Initial Guess
        beta, r, Ay = x   #sideslip (rad), yaw rate (rad/s), lateral accel (m/s^2)

        #BLOCK for it in range(max_iter):
        #Slip Angles (Rads)
        Vy = Vx * beta
        alpha_FL = (Vy+r*a)/(Vx-r*(vp.FTrackwidth_mm/1000/2)) - delta #Different delta if ackerman
        alpha_FR = (Vy+r*a)/(Vx+r*(vp.FTrackwidth_mm/1000/2)) - delta
        alpha_RL = (Vy-r*b)/(Vx-r*(vp.RTrackwidth_mm/1000/2))
        alpha_RR = (Vy-r*b)/(Vx+r*(vp.RTrackwidth_mm/1000/2))

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

        #Will add Lateral Load Transfer here later

        FZ_FL = FZ_FL + DF_FL
        FZ_FR = FZ_FR + DF_FR
        FZ_RL = FZ_RL + DF_RL
        FZ_RR = FZ_RR + DF_RR

        #Find FY
        gamma_FL = vp.Camber_By_Travel_deg(0)
        gamma_FR = vp.Camber_By_Travel_deg(0)
        gamma_RL = vp.Camber_By_Travel_deg(0)
        gamma_RR = vp.Camber_By_Travel_deg(0)

        FY_FL = tire.FY(alpha_FL*RAD2DEG, FZ_FL*N2LBF, gamma_FL)*LBF2N
        FY_FR = tire.FY(alpha_FR*RAD2DEG, FZ_FR*N2LBF, gamma_FR)*LBF2N
        FY_RL = tire.FY(alpha_RL*RAD2DEG, FZ_RL*N2LBF, gamma_RL)*LBF2N
        FY_RR = tire.FY(alpha_RR*RAD2DEG, FZ_RR*N2LBF, gamma_RR)*LBF2N

        MZ_FL = tire.MZ(alpha_FL*RAD2DEG, FZ_FL*N2LBF, gamma_FL)*FTLB2NM
        MZ_FR = tire.MZ(alpha_FR*RAD2DEG, FZ_FR*N2LBF, gamma_FR)*FTLB2NM
        MZ_RL = tire.MZ(alpha_RL*RAD2DEG, FZ_RL*N2LBF, gamma_RL)*FTLB2NM
        MZ_RR = tire.MZ(alpha_RR*RAD2DEG, FZ_RR*N2LBF, gamma_RR)*FTLB2NM

        #Force-Moment Balance
        FYF = FY_FL + FY_FR
        FYR = FY_RL + FY_RR

        MZ_Total = MZ_FL + MZ_FR + MZ_RL + MZ_RR

        #Make Residuals
        R1 = m*Ay - (FYF*np.cos(delta)+FYR)
        R2 = a*FYF*np.cos(delta) - b*FYR + MZ_Total
        R3 = Ay - Vx*r

        #BLOCK Converge plz
        #d_beta = relax*(beta_new-beta)
        #d_r = relax*(r_new-r)
        #d_Ay = relax*(Ay_new-Ay)
        #beta += d_beta
        #r += d_r
        #Ay += d_Ay
        #if abs(d_beta) < tol_beta and abs(d_r) < tol_r and abs(d_Ay) < tol_Ay:
        #    break
        #else:
        #    raise RuntimeError("Fella did not converge")

        return [R1, R2, R3]
    result = ls(residual, np.array([0.0, 0.0, 0.0]))
    if result.cost > 1e-6:
        raise RuntimeError(
            f"Solver failed to converge: cost={result.cost:.2e}"
        )

    beta, r, Ay = result.x
    return beta, r, Ay











