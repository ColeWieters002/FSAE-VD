import os
import numpy as np
from scipy.interpolate import interp1d as curve


_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_TIRE_DIR = os.path.join(_THIS_DIR, "Tires")


# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
LBF2N = 4.4482216
N2LBF = 1.0 / LBF2N
FTLB2NM = 1.3558
NM2FTLB = 1.0 / FTLB2NM
FT2M = 0.3048
M2FT = 1.0 / FT2M
IN2M = 0.0254
M2IN = 1.0 / IN2M
RAD2DEG = 180.0 / np.pi
DEG2RAD = np.pi / 180.0
MPS2MPH = 2.23694
MPH2MPS = 1.0 / MPS2MPH


# ---------------------------------------------------------------------------
# Actuation / stiffness
# ---------------------------------------------------------------------------
FrontRollStiffness = 25000.0 * FTLB2NM   # N*m/rad
RearRollStiffness = 23000.0 * FTLB2NM    # N*m/rad
FrontHeaveStiffness = 400.0 * LBF2N / IN2M  # N/m, axle heave rate
RearHeaveStiffness = 375.0 * LBF2N / IN2M   # N/m, axle heave rate


# ---------------------------------------------------------------------------
# Mass / hardpoints
# ---------------------------------------------------------------------------
Vehicle_kg = 176.0
Driver_kg = 76.0
TotalMass_kg = Vehicle_kg + Driver_kg
UnsprungMass_kg = 42.0
SprungMass_kg = TotalMass_kg - UnsprungMass_kg
WeightDist = 0.48  # front fraction
CG_mm = 11.25 * IN2M * 1000.0
YawInertia = 92.0  # kg*m^2; unused by the quasi-static YMD

Wheelbase_mm = 60.5 * IN2M * 1000.0
FTrackwidth_mm = 50.0 * IN2M * 1000.0
RTrackwidth_mm = 49.0 * IN2M * 1000.0


# ---------------------------------------------------------------------------
# Camber geometry
#
# Convention used below:
#   * A negative wheel-local camber angle is negative camber.
#   * The right wheel's vehicle-coordinate camber is the sign-mirror of its
#     wheel-local value, because the supplied TIR is a LEFT-tire model.
#   * Positive body roll is the YMD convention used here: right side outside.
# ---------------------------------------------------------------------------

# Keep the existing design's nominal static camber target.  These replace the
# old call Camber_By_Travel_deg(-2), which accidentally meant -2 mm rather
# than a clear static-camber setting.
FrontStaticCamber_deg = -2.0
RearStaticCamber_deg = -2.0

# Measured geometry from the supplied camber-gain slide.
FrontHeaveCamberGain_deg_per_mm = -0.0425
RearHeaveCamberGain_deg_per_mm = -0.0471
FrontRollCamberGain_deg_per_deg = -0.528
RearRollCamberGain_deg_per_deg = -0.488
FrontSteerCamberGain_deg_per_deg = -0.0275
RearSteerCamberGain_deg_per_deg = 0.0

# The YMD has no absolute heave state.  Its static camber is therefore defined
# at static ride height, and only aero-induced compression is added as heave.
IncludeAeroHeaveCamber = True


def CamberAtCorner_deg(axle, side, roll_deg, steer_deg=0.0, heave_mm=0.0):
    """Return vehicle-coordinate camber for one wheel.

    Parameters
    ----------
    axle : "front" or "rear"
    side : "left" or "right"
    roll_deg : body roll, positive when the right tyre is outside
    steer_deg : road-wheel steer; used only at the front
    heave_mm : positive bump/compression from static ride height

    The available steering result is a scalar inside-tire gain, not separate
    left/right fits.  Applying it against |steer| to both front tyre-local
    cambers keeps a left/right turn symmetric until wheel-specific data is
    available.
    """
    axle = axle.lower()
    side = side.lower()
    if axle not in ("front", "rear"):
        raise ValueError("axle must be 'front' or 'rear'")
    if side not in ("left", "right"):
        raise ValueError("side must be 'left' or 'right'")

    if axle == "front":
        static_camber = FrontStaticCamber_deg
        heave_gain = FrontHeaveCamberGain_deg_per_mm
        roll_gain = FrontRollCamberGain_deg_per_deg
        steer_gain = FrontSteerCamberGain_deg_per_deg
    else:
        static_camber = RearStaticCamber_deg
        heave_gain = RearHeaveCamberGain_deg_per_mm
        roll_gain = RearRollCamberGain_deg_per_deg
        steer_gain = RearSteerCamberGain_deg_per_deg

    # Positive roll -> right is outside -> its local camber becomes more
    # negative because roll_gain is negative; left becomes less negative.
    roll_term = roll_gain * roll_deg if side == "right" else -roll_gain * roll_deg
    local_camber = (
        static_camber
        + heave_gain * heave_mm
        + roll_term
        + steer_gain * abs(steer_deg)
    )

    # Convert local model convention back to vehicle-coordinate camber.
    return local_camber if side == "left" else -local_camber


# Retained only for older scripts that still call this function.  The repaired
# YMD uses CamberAtCorner_deg above rather than this generic placeholder curve.
CamberBounds = [0.0, -2.0, -3.0]
_CamberCurve = curve([-25.4, 0.0, 25.4], CamberBounds, kind="quadratic", fill_value="extrapolate")


def Camber_By_Travel_deg(travel_mm, side):
    camber = float(_CamberCurve(travel_mm))
    if side.lower() == "left":
        return camber
    if side.lower() == "right":
        return -camber
    raise ValueError("side must be 'left' or 'right'")


# ---------------------------------------------------------------------------
# Suspension / steering geometry
# ---------------------------------------------------------------------------
FrontRollCenter_mm = -0.02 * IN2M * 1000.0
RearRollCenter_mm = 2.474 * IN2M * 1000.0
Ackerman = 0.0
FrontStaticToe_deg = 0.0
RearStaticToe_deg = 0.0
FrontCaster_deg = 8.23
RearCaster_deg = 2.81
FrontKPI_deg = 0.0
RearKPI_deg = 14.9


# ---------------------------------------------------------------------------
# Environment / aero / tyre
# ---------------------------------------------------------------------------
Gravity = 9.8
AirDensity = 1.225
A = 1.0
CL = -3.75
CD = 1.4
AeroBalance = 0.4  # front fraction
TirePressure_bar = 0.827  # 12 psi
RollResistance = 0.4

Hoosier_16x75_10_R20 = os.path.join(_TIRE_DIR, "Hoosier_16x75_10_R20.tir")
Hoosier_18x75_10_R20 = os.path.join(_TIRE_DIR, "Hoosier_18x75_10_R20.tir")
Hoosier_16x75_10_LC0 = os.path.join(_TIRE_DIR, "Hoosier_16x75_10_LC0.tir")
TireModel = Hoosier_16x75_10_R20
