import os
import csv
import numpy as np
import matplotlib.pyplot as plt

import VehicleParameters as vp
from VehicleParameters import LBF2N, N2LBF, FTLB2NM
from YMDSim import Tire


# ============================================================
# USER SETTINGS
# ============================================================

TIRE_FILE = vp.TireModel
TIRE_PRESSURE_BAR = vp.TirePressure_bar

SHOW_PLOTS = True
SAVE_PLOTS = True
SAVE_CSV = True
OUTPUT_FOLDER = "TireAnalysis_Output"

# All loads entered below are NEWTONS.
# Slip angle and camber are DEGREES.
# Slip ratio is a decimal: 0.10 = 10% slip.

SWEEPS = [
    {
        "name": "Slip Angle vs Lateral Force - Camber Sweep",
        "x": "SA",
        "y": "FY",
        "sweep": "IA",
        "x_values": np.linspace(-15.0, 15.0, 301),
        "sweep_values": [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0],
        "fixed": {
            "FZ": 700.0,
            "SR": 0.0,
        },
    },
    {
        "name": "Slip Angle vs Lateral Force - Load Sweep",
        "x": "SA",
        "y": "FY",
        "sweep": "FZ",
        "x_values": np.linspace(-15.0, 15.0, 301),
        "sweep_values": [300.0, 500.0, 700.0, 900.0, 1100.0],
        "fixed": {
            "IA": 0.0,
            "SR": 0.0,
        },
    },
    {
        "name": "Slip Angle vs Aligning Moment - Camber Sweep",
        "x": "SA",
        "y": "MZ",
        "sweep": "IA",
        "x_values": np.linspace(-15.0, 15.0, 301),
        "sweep_values": [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0],
        "fixed": {
            "FZ": 700.0,
            "SR": 0.0,
        },
    },
    {
        "name": "Slip Ratio vs Longitudinal Force - Load Sweep",
        "x": "SR",
        "y": "FX",
        "sweep": "FZ",
        "x_values": np.linspace(-0.30, 0.30, 301),
        "sweep_values": [300.0, 500.0, 700.0, 900.0, 1100.0],
        "fixed": {
            "IA": 0.0,
            "SA": 0.0,
        },
    },
    {
        "name": "Slip Angle vs Lateral Friction Coefficient - Load Sweep",
        "x": "SA",
        "y": "MUY",
        "sweep": "FZ",
        "x_values": np.linspace(-15.0, 15.0, 301),
        "sweep_values": [300.0, 500.0, 700.0, 900.0, 1100.0],
        "fixed": {
            "IA": 0.0,
            "SR": 0.0,
        },
    },
]


# ============================================================
# INTERNAL DEFINITIONS
# ============================================================

DEFAULT_INPUTS = {
    "SA": 0.0,     # Slip angle [deg]
    "SR": 0.0,     # Slip ratio [-]
    "IA": 0.0,     # Inclination / camber angle [deg]
    "FZ": 700.0,   # Vertical load [N]
}

AXIS_LABELS = {
    "SA": "Slip Angle (deg)",
    "SR": "Slip Ratio",
    "IA": "Camber Angle (deg)",
    "FZ": "Vertical Load (N)",
    "FY": "Lateral Force (N)",
    "FX": "Longitudinal Force (N)",
    "MZ": "Self-Aligning Moment (N*m)",
    "MUY": "Lateral Friction Coefficient |Fy| / Fz",
    "MUX": "Longitudinal Friction Coefficient |Fx| / Fz",
    "TRAIL": "Effective Pneumatic Trail (mm)",
}


# ============================================================
# TIRE MODEL EVALUATION
# ============================================================

def evaluate_tire(tire, output, SA=0.0, SR=0.0, IA=0.0, FZ=700.0):
    """
    User-facing units:
        SA = slip angle [deg]
        SR = slip ratio [-]
        IA = camber/inclination angle [deg]
        FZ = vertical load [N]

    Outputs:
        FY = lateral force [N]
        FX = longitudinal force [N]
        MZ = aligning moment [N*m]
        MUY = |FY| / FZ [-]
        MUX = |FX| / FZ [-]
        TRAIL = -MZ/FY [mm]
    """

    FZ_N = np.asarray(FZ, dtype=float)
    FZ_lbf = FZ_N * N2LBF

    if output == "FY":
        return np.asarray(tire.FY(SA, FZ_lbf, IA), dtype=float) * LBF2N

    if output == "FX":
        return np.asarray(tire.FX(SR, FZ_lbf, IA), dtype=float) * LBF2N

    if output == "MZ":
        return np.asarray(tire.MZ(SA, FZ_lbf, IA), dtype=float) * FTLB2NM

    if output == "MUY":
        FY_N = evaluate_tire(tire, "FY", SA=SA, SR=SR, IA=IA, FZ=FZ_N)
        return np.divide(
            np.abs(FY_N),
            FZ_N,
            out=np.full(np.broadcast(FY_N, FZ_N).shape, np.nan, dtype=float),
            where=np.asarray(FZ_N) != 0.0,
        )

    if output == "MUX":
        FX_N = evaluate_tire(tire, "FX", SA=SA, SR=SR, IA=IA, FZ=FZ_N)
        return np.divide(
            np.abs(FX_N),
            FZ_N,
            out=np.full(np.broadcast(FX_N, FZ_N).shape, np.nan, dtype=float),
            where=np.asarray(FZ_N) != 0.0,
        )

    if output == "TRAIL":
        FY_N = evaluate_tire(tire, "FY", SA=SA, SR=SR, IA=IA, FZ=FZ_N)
        MZ_NM = evaluate_tire(tire, "MZ", SA=SA, SR=SR, IA=IA, FZ=FZ_N)

        trail_m = np.divide(
            -MZ_NM,
            FY_N,
            out=np.full(np.broadcast(MZ_NM, FY_N).shape, np.nan, dtype=float),
            where=np.abs(FY_N) > 1.0,
        )

        return trail_m * 1000.0

    raise ValueError(f"Unknown output '{output}'")


# ============================================================
# SWEEP TOOLS
# ============================================================

def build_inputs(x_name, x_values, sweep_name, sweep_value, fixed):
    inputs = DEFAULT_INPUTS.copy()
    inputs.update(fixed)

    inputs[x_name] = x_values
    inputs[sweep_name] = sweep_value

    return inputs


def format_value(variable, value):
    if variable in ["SA", "IA"]:
        return f"{value:.1f} deg"

    if variable == "FZ":
        return f"{value:.0f} N"

    if variable == "SR":
        return f"{100.0 * value:.1f}%"

    return f"{value}"


def safe_filename(name):
    keep = []

    for character in name:
        if character.isalnum() or character in ["-", "_"]:
            keep.append(character)
        elif character in [" ", "/"]:
            keep.append("_")

    return "".join(keep)


def print_curve_summary(x_name, y_name, x_values, y_values, inputs, sweep_name, sweep_value):
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)

    finite = np.isfinite(x_values) & np.isfinite(y_values)

    if not np.any(finite):
        return

    x = x_values[finite]
    y = y_values[finite]

    index_peak = np.argmax(np.abs(y))
    x_peak = x[index_peak]
    y_peak = y[index_peak]

    label = f"{sweep_name}={format_value(sweep_name, sweep_value)}"

    if y_name == "FY":
        FZ = np.asarray(inputs["FZ"], dtype=float)
        FZ_peak = FZ if FZ.ndim == 0 else FZ[finite][index_peak]
        mu_peak = abs(y_peak) / FZ_peak if FZ_peak != 0.0 else np.nan

        message = (
            f"  {label:<22} "
            f"Peak |Fy|={abs(y_peak):8.1f} N at {x_name}={x_peak:7.3f}, "
            f"mu_y={mu_peak:5.3f}"
        )

        if x_name == "SA" and len(x) >= 3:
            dF_dSA = np.gradient(y, x)
            zero_index = np.argmin(np.abs(x))
            cornering_stiffness = abs(dF_dSA[zero_index]) * 180.0 / np.pi
            message += f", |C_alpha|={cornering_stiffness:9.1f} N/rad"

        print(message)
        return

    if y_name == "FX":
        FZ = np.asarray(inputs["FZ"], dtype=float)
        FZ_peak = FZ if FZ.ndim == 0 else FZ[finite][index_peak]
        mu_peak = abs(y_peak) / FZ_peak if FZ_peak != 0.0 else np.nan

        print(
            f"  {label:<22} "
            f"Peak |Fx|={abs(y_peak):8.1f} N at {x_name}={x_peak:7.3f}, "
            f"mu_x={mu_peak:5.3f}"
        )
        return

    if y_name == "MZ":
        print(
            f"  {label:<22} "
            f"Peak |Mz|={abs(y_peak):8.2f} N*m at {x_name}={x_peak:7.3f}"
        )
        return

    if y_name in ["MUY", "MUX"]:
        print(
            f"  {label:<22} "
            f"Peak {y_name}={np.nanmax(y):6.3f}"
        )


def save_csv(plot_config, x_values, curves):
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    filename = safe_filename(plot_config["name"]) + ".csv"
    filepath = os.path.join(OUTPUT_FOLDER, filename)

    headers = [AXIS_LABELS[plot_config["x"]]]
    headers += [curve["label"] for curve in curves]

    with open(filepath, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headers)

        for i in range(len(x_values)):
            row = [x_values[i]]
            row += [curve["y"][i] for curve in curves]
            writer.writerow(row)


def plot_sweep(tire, plot_config):
    x_name = plot_config["x"]
    y_name = plot_config["y"]
    sweep_name = plot_config["sweep"]

    if x_name == sweep_name:
        raise ValueError("The x variable and sweep variable cannot be the same")

    x_values = np.asarray(plot_config["x_values"], dtype=float)
    sweep_values = plot_config["sweep_values"]
    fixed = plot_config.get("fixed", {})

    plt.figure(figsize=(10, 7))

    curves = []

    print()
    print("=" * 90)
    print(plot_config["name"])
    print("=" * 90)

    for sweep_value in sweep_values:
        inputs = build_inputs(
            x_name=x_name,
            x_values=x_values,
            sweep_name=sweep_name,
            sweep_value=sweep_value,
            fixed=fixed,
        )

        y_values = evaluate_tire(tire, y_name, **inputs)

        # If the tire function returned a scalar, make it a full-length curve.
        if np.asarray(y_values).ndim == 0:
            y_values = np.full_like(x_values, float(y_values), dtype=float)

        label = f"{sweep_name} = {format_value(sweep_name, sweep_value)}"

        plt.plot(x_values, y_values, label=label)

        curves.append({
            "label": label,
            "y": np.asarray(y_values, dtype=float),
        })

        print_curve_summary(
            x_name,
            y_name,
            x_values,
            y_values,
            inputs,
            sweep_name,
            sweep_value,
        )

    plt.title(plot_config["name"])
    plt.xlabel(AXIS_LABELS[x_name])
    plt.ylabel(AXIS_LABELS[y_name])
    plt.axhline(0.0, linewidth=1.0)

    if np.nanmin(x_values) <= 0.0 <= np.nanmax(x_values):
        plt.axvline(0.0, linewidth=1.0)

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    if SAVE_PLOTS:
        os.makedirs(OUTPUT_FOLDER, exist_ok=True)
        filename = safe_filename(plot_config["name"]) + ".png"
        filepath = os.path.join(OUTPUT_FOLDER, filename)
        plt.savefig(filepath, dpi=200)

    if SAVE_CSV:
        save_csv(plot_config, x_values, curves)

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


# ============================================================
# MAIN
# ============================================================

def main():
    tire = Tire(TIRE_FILE, TIRE_PRESSURE_BAR)

    print()
    print("=" * 90)
    print("TIRE CHARACTERIZATION")
    print("=" * 90)
    print(f"Tire file       : {TIRE_FILE}")
    print(f"Pressure stored : {TIRE_PRESSURE_BAR:.3f} bar")
    print(f"Nominal load    : {tire.fz0:.2f} lbf = {tire.fz0 * LBF2N:.2f} N")
    print()
    print("NOTE: Your current MF5.2 FY/FX/MZ functions do not use pressure as an input.")
    print("Changing TIRE_PRESSURE_BAR will therefore not change these curves yet.")

    for plot_config in SWEEPS:
        plot_sweep(tire, plot_config)


if __name__ == "__main__":
    main()
