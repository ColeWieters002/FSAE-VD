"""Quasi-static dual-track YMD force solver.

This is intentionally force-equivalent to the supplied solver.  The only
functional cleanup is a safe LLTD diagnostic at Ay = 0, where total load
transfer is zero and LLTD is undefined.
"""

import numpy as np
from scipy.optimize import least_squares as ls

import TireFunctions as MF
import VehicleParameters as vp
from VehicleParameters import DEG2RAD, RAD2DEG, LBF2N, N2LBF, FT2M


class Tire:
    def __init__(self, tir_path, pressure_bar=None):
        sections = self._parse_tir(tir_path)
        self.p_fy = np.array(
            [sections["LATERAL_COEFFICIENTS"][name] for name in MF.FY_Params],
            dtype=float,
        )
        self.q_mz = np.array(
            [sections["ALIGNING_COEFFICIENTS"][name] for name in MF.MZ_Params],
            dtype=float,
        )
        self.p_fx = np.array(
            [sections["LONGITUDINAL_COEFFICIENTS"][name] for name in MF.FX_Params],
            dtype=float,
        )
        self.fz0 = sections["VERTICAL"]["FNOMIN"] / LBF2N
        self.r0 = sections["DIMENSION"]["UNLOADED_RADIUS"] / FT2M
        self.pressure_bar = pressure_bar
        self.path = tir_path

    def FY(self, slip_angle_deg, fz_lbf, camber_deg):
        return MF.FY((slip_angle_deg, fz_lbf, camber_deg), self.p_fy, self.fz0)

    def MZ(self, slip_angle_deg, fz_lbf, camber_deg):
        return MF.MZ(
            (slip_angle_deg, fz_lbf, camber_deg),
            self.q_mz,
            self.p_fy,
            self.fz0,
            self.r0,
        )

    def FX(self, slip_ratio, fz_lbf, camber_deg):
        return MF.FX((slip_ratio, fz_lbf, camber_deg), self.p_fx, self.fz0)

    @staticmethod
    def _parse_tir(path):
        sections = {}
        current_section = None
        with open(path, "r", encoding="utf-8") as tire_file:
            for raw_line in tire_file:
                line = raw_line.split("!", 1)[0].strip()
                if not line:
                    continue
                if line.startswith("[") and line.endswith("]"):
                    current_section = line[1:-1].strip().upper()
                    sections[current_section] = {}
                    continue
                if current_section is None or "=" not in line:
                    continue
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip().strip("'").strip('"')
                try:
                    value = float(value)
                except ValueError:
                    pass
                sections[current_section][key] = value
        return sections


def solve(vx_mps, beta_rad, delta_rad, vehicle=vp, tire=None, debug=False):
    """Solve lateral-force and roll equilibrium for one beta/delta point.

    The yaw moment is deliberately *not* constrained to zero; that is what
    makes this a yaw-moment-diagram point rather than a steady-state trim
    solution.
    """
    if tire is None:
        tire = Tire(vehicle.TireModel, vehicle.TirePressure_bar)

    wheelbase = vehicle.Wheelbase_mm / 1000.0
    a = wheelbase * (1.0 - vehicle.WeightDist)
    b = wheelbase * vehicle.WeightDist
    mass = vehicle.TotalMass_kg
    sprung_mass = vehicle.SprungMass_kg
    gravity = vehicle.Gravity
    front_track = vehicle.FTrackwidth_mm / 1000.0
    rear_track = vehicle.RTrackwidth_mm / 1000.0
    front_rc = vehicle.FrontRollCenter_mm / 1000.0
    rear_rc = vehicle.RearRollCenter_mm / 1000.0
    roll_axis_height = (front_rc * b + rear_rc * a) / wheelbase
    cg_height = vehicle.CG_mm / 1000.0
    total_roll_stiffness = vehicle.FrontRollStiffness + vehicle.RearRollStiffness

    def calculate_state(ay_mps2, roll_rad):
        # This is the r = 0 instantaneous/quasi-static YMD formulation.
        yaw_rate = 0.0
        vy_mps = vx_mps * np.tan(beta_rad)

        alpha_fl = delta_rad - np.arctan2(vy_mps + yaw_rate * a, vx_mps - yaw_rate * front_track / 2.0)
        alpha_fr = delta_rad - np.arctan2(vy_mps + yaw_rate * a, vx_mps + yaw_rate * front_track / 2.0)
        alpha_rl = -np.arctan2(vy_mps - yaw_rate * b, vx_mps - yaw_rate * rear_track / 2.0)
        alpha_rr = -np.arctan2(vy_mps - yaw_rate * b, vx_mps + yaw_rate * rear_track / 2.0)

        # Static corner weights.
        fz_fl = mass * gravity * vehicle.WeightDist / 2.0
        fz_fr = mass * gravity * vehicle.WeightDist / 2.0
        fz_rl = mass * gravity * (1.0 - vehicle.WeightDist) / 2.0
        fz_rr = mass * gravity * (1.0 - vehicle.WeightDist) / 2.0

        # CL is negative for downforce in the supplied parameter convention.
        aero_force = 0.5 * vehicle.AirDensity * vx_mps**2 * vehicle.CL * vehicle.A
        df_fl = -aero_force * vehicle.AeroBalance * 0.5
        df_fr = -aero_force * vehicle.AeroBalance * 0.5
        df_rl = -aero_force * (1.0 - vehicle.AeroBalance) * 0.5
        df_rr = -aero_force * (1.0 - vehicle.AeroBalance) * 0.5

        # Static camber is defined at static ride height.  The YMD has no
        # absolute heave state, so only aero's added vertical load produces a
        # symmetric heave displacement about that condition.
        if getattr(vehicle, "IncludeAeroHeaveCamber", True):
            front_aero_downforce_n = -aero_force * vehicle.AeroBalance
            rear_aero_downforce_n = -aero_force * (1.0 - vehicle.AeroBalance)
            front_heave_mm = 1000.0 * front_aero_downforce_n / vehicle.FrontHeaveStiffness
            rear_heave_mm = 1000.0 * rear_aero_downforce_n / vehicle.RearHeaveStiffness
        else:
            front_heave_mm = 0.0
            rear_heave_mm = 0.0

        # Quasi-static geometric + elastic lateral load transfer.
        sprung_lateral_force = sprung_mass * ay_mps2
        front_lateral_share = sprung_lateral_force * b / wheelbase
        rear_lateral_share = sprung_lateral_force * a / wheelbase
        front_geo_moment = front_lateral_share * front_rc
        rear_geo_moment = rear_lateral_share * rear_rc
        front_elastic_moment = vehicle.FrontRollStiffness * roll_rad
        rear_elastic_moment = vehicle.RearRollStiffness * roll_rad
        front_roll_moment = front_geo_moment + front_elastic_moment
        rear_roll_moment = rear_geo_moment + rear_elastic_moment
        front_load_transfer = front_roll_moment / front_track
        rear_load_transfer = rear_roll_moment / rear_track

        df_fl -= front_load_transfer
        df_fr += front_load_transfer
        df_rl -= rear_load_transfer
        df_rr += rear_load_transfer

        fz_fl_raw = fz_fl + df_fl
        fz_fr_raw = fz_fr + df_fr
        fz_rl_raw = fz_rl + df_rl
        fz_rr_raw = fz_rr + df_rr

        # Preserve axle vertical load if an inside wheel lifts.
        fz_fl, fz_fr = fz_fl_raw, fz_fr_raw
        if fz_fl < 0.0:
            fz_fr += fz_fl
            fz_fl = 0.0
        elif fz_fr < 0.0:
            fz_fl += fz_fr
            fz_fr = 0.0

        fz_rl, fz_rr = fz_rl_raw, fz_rr_raw
        if fz_rl < 0.0:
            fz_rr += fz_rl
            fz_rl = 0.0
        elif fz_rr < 0.0:
            fz_rl += fz_rr
            fz_rr = 0.0

        fz_fl = max(fz_fl, 0.0)
        fz_fr = max(fz_fr, 0.0)
        fz_rl = max(fz_rl, 0.0)
        fz_rr = max(fz_rr, 0.0)

        # Measured heave, roll, and steering camber gains.  The vehicle
        # parameter helper returns vehicle-coordinate camber; right tyre
        # inputs are mirrored below for the left-tyre TIR convention.
        roll_deg = roll_rad * RAD2DEG
        steer_deg = delta_rad * RAD2DEG
        gamma_fl = vehicle.CamberAtCorner_deg("front", "left", roll_deg, steer_deg, front_heave_mm)
        gamma_fr = vehicle.CamberAtCorner_deg("front", "right", roll_deg, steer_deg, front_heave_mm)
        gamma_rl = vehicle.CamberAtCorner_deg("rear", "left", roll_deg, 0.0, rear_heave_mm)
        gamma_rr = vehicle.CamberAtCorner_deg("rear", "right", roll_deg, 0.0, rear_heave_mm)

        # Equivalent wheel travel relative to static ride height.  Positive is
        # bump/compression; negative is rebound/extension.  This is a rigid
        # body roll approximation at the wheel centerline, not spring travel,
        # so a motion-ratio conversion is still required for damper stroke.
        # Positive roll in this YMD convention loads the right-hand tyres.
        front_roll_travel_mm = 1000.0 * (front_track / 2.0) * np.tan(roll_rad)
        rear_roll_travel_mm = 1000.0 * (rear_track / 2.0) * np.tan(roll_rad)
        travel_fl_mm = front_heave_mm - front_roll_travel_mm
        travel_fr_mm = front_heave_mm + front_roll_travel_mm
        travel_rl_mm = rear_heave_mm - rear_roll_travel_mm
        travel_rr_mm = rear_heave_mm + rear_roll_travel_mm

        # Mirror the left-tyre TIR model for the right side.
        fy_fl = tire.FY(alpha_fl * RAD2DEG, fz_fl * N2LBF, gamma_fl) * LBF2N
        fy_fr = -tire.FY(-alpha_fr * RAD2DEG, fz_fr * N2LBF, -gamma_fr) * LBF2N
        fy_rl = tire.FY(alpha_rl * RAD2DEG, fz_rl * N2LBF, gamma_rl) * LBF2N
        fy_rr = -tire.FY(-alpha_rr * RAD2DEG, fz_rr * N2LBF, -gamma_rr) * LBF2N

        # Rotate front lateral forces from wheel to body axes.
        fx_fl_body = -fy_fl * np.sin(delta_rad)
        fy_fl_body = fy_fl * np.cos(delta_rad)
        fx_fr_body = -fy_fr * np.sin(delta_rad)
        fy_fr_body = fy_fr * np.cos(delta_rad)
        fx_rl_body = 0.0
        fy_rl_body = fy_rl
        fx_rr_body = 0.0
        fy_rr_body = fy_rr

        moment_fl = a * fy_fl_body - (front_track / 2.0) * fx_fl_body
        moment_fr = a * fy_fr_body - (-front_track / 2.0) * fx_fr_body
        moment_rl = -b * fy_rl_body
        moment_rr = -b * fy_rr_body
        force_moment_total = moment_fl + moment_fr + moment_rl + moment_rr
        fy_total = fy_fl_body + fy_fr_body + fy_rl_body + fy_rr_body

        total_load_transfer = front_load_transfer + rear_load_transfer
        if abs(total_load_transfer) < 1e-12:
            tlltd_front = np.nan
            tlltd_rear = np.nan
        else:
            tlltd_front = front_load_transfer / total_load_transfer
            tlltd_rear = 1.0 - tlltd_front

        data = {
            "Ay": ay_mps2,
            "Ay_g": ay_mps2 / gravity,
            "phi_rad": roll_rad,
            "phi_deg": roll_rad * RAD2DEG,
            "FrontLoadTransfer_N": front_load_transfer,
            "RearLoadTransfer_N": rear_load_transfer,
            "TLLTD_Front": tlltd_front,
            "TLLTD_Rear": tlltd_rear,
            "RollAxisHeight_m": roll_axis_height,
            "FrontGeoMoment_Nm": front_geo_moment,
            "RearGeoMoment_Nm": rear_geo_moment,
            "FrontElasticMoment_Nm": front_elastic_moment,
            "RearElasticMoment_Nm": rear_elastic_moment,
            "FrontRollMoment_Nm": front_roll_moment,
            "RearRollMoment_Nm": rear_roll_moment,
            "FrontAeroHeave_mm": front_heave_mm,
            "RearAeroHeave_mm": rear_heave_mm,
            "FrontRollTravel_mm": front_roll_travel_mm,
            "RearRollTravel_mm": rear_roll_travel_mm,
            "Travel_FL_mm": travel_fl_mm,
            "Travel_FR_mm": travel_fr_mm,
            "Travel_RL_mm": travel_rl_mm,
            "Travel_RR_mm": travel_rr_mm,
            "gamma_FL_deg": gamma_fl,
            "gamma_FR_deg": gamma_fr,
            "gamma_RL_deg": gamma_rl,
            "gamma_RR_deg": gamma_rr,
            "FZ_FL_N": fz_fl,
            "FZ_FR_N": fz_fr,
            "FZ_RL_N": fz_rl,
            "FZ_RR_N": fz_rr,
            "FZ_FL_raw_N": fz_fl_raw,
            "FZ_FR_raw_N": fz_fr_raw,
            "FZ_RL_raw_N": fz_rl_raw,
            "FZ_RR_raw_N": fz_rr_raw,
            "FY_FL_N": fy_fl,
            "FY_FR_N": fy_fr,
            "FY_RL_N": fy_rl,
            "FY_RR_N": fy_rr,
            "FY_Total_N": fy_total,
            "alpha_FL_deg": alpha_fl * RAD2DEG,
            "alpha_FR_deg": alpha_fr * RAD2DEG,
            "alpha_RL_deg": alpha_rl * RAD2DEG,
            "alpha_RR_deg": alpha_rr * RAD2DEG,
            "ForceMoment_FL_Nm": moment_fl,
            "ForceMoment_FR_Nm": moment_fr,
            "ForceMoment_RL_Nm": moment_rl,
            "ForceMoment_RR_Nm": moment_rr,
            "ForceMoment_Total_Nm": force_moment_total,
            "Mz_Vehicle_Nm": force_moment_total,
        }

        if debug:
            _append_debug(vx_mps, beta_rad, delta_rad, data)

        return fy_total, force_moment_total, data

    def residual(state):
        ay_mps2, roll_rad = state
        fy_total, _, _ = calculate_state(ay_mps2, roll_rad)
        lateral_force_residual = ay_mps2 - fy_total / mass
        roll_residual = total_roll_stiffness * roll_rad - sprung_mass * ay_mps2 * (cg_height - roll_axis_height)
        return [lateral_force_residual, roll_residual]

    result = ls(residual, np.array([0.0, 0.0]))
    if result.cost > 1e-6:
        print(f"warning: residual not driven to zero (cost={result.cost:.2e})")

    ay_mps2, roll_rad = result.x
    fy_total, yaw_moment_nm, data = calculate_state(ay_mps2, roll_rad)
    return ay_mps2, yaw_moment_nm, roll_rad, data


def _append_debug(vx_mps, beta_rad, delta_rad, data):
    with open("YMD_Debug.txt", "a", encoding="utf-8") as debug_file:
        debug_file.write("\n====================================\n")
        debug_file.write(f"Vx = {vx_mps}\n")
        debug_file.write(f"beta = {beta_rad * RAD2DEG} deg\n")
        debug_file.write(f"delta = {delta_rad * RAD2DEG} deg\n")
        for heading, keys in (
            ("Slip Angles (deg)", ("alpha_FL_deg", "alpha_FR_deg", "alpha_RL_deg", "alpha_RR_deg")),
            ("FZ (N)", ("FZ_FL_N", "FZ_FR_N", "FZ_RL_N", "FZ_RR_N")),
            ("FY (N)", ("FY_FL_N", "FY_FR_N", "FY_RL_N", "FY_RR_N")),
            ("Wheel Travel (mm; +compression, -extension)", ("Travel_FL_mm", "Travel_FR_mm", "Travel_RL_mm", "Travel_RR_mm")),
        ):
            debug_file.write(f"\n{heading}\n")
            for key in keys:
                debug_file.write(f"{key} = {data[key]}\n")
        debug_file.write(f"FY Total = {data['FY_Total_N']}\n")
        debug_file.write(f"TOTAL Mz = {data['Mz_Vehicle_Nm']}\n")
