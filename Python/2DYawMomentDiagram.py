"""2-D yaw-moment diagram with a physically usable operating-envelope plot.

The raw YMD is still calculated at every beta/delta point.  The plotted map
keeps only the continuous branch that contains beta = 0 and has monotonically
decreasing Ay as beta increases.  Once a tyre-force combination goes
post-peak, Ay reverses direction and the raw parametric line curls back on
itself; that branch is intentionally excluded from the operating plot.
"""

import numpy as np
import matplotlib.pyplot as plt

import VehicleParameters as vp
from VehicleParameters import DEG2RAD, RAD2DEG
from YMDSim import Tire, solve


# ---------------------------------------------------------------------------
# User inputs
# ---------------------------------------------------------------------------
VX_MPS = 45 / vp.MPS2MPH
SWEEP_LIMIT_DEG = 13.0
GRID_STEP_DEG = 1

# This is an explicit operating limit, not a claim about the tire model's
# validated range.  Set to None to use only the monotonic-branch criterion.
MAX_TIRE_SLIP_DEG = 12.0

DEBUG = False
PLOT_WHOLE_DEGREE_LINES_ONLY = True


def is_whole_degree(value):
    return abs(value - round(value)) < 1e-9


def monotonic_branch_mask(ay_grid, beta_values):
    """Keep the no-turnback portion of each constant-steer line.

    On the operating branch, increasing beta decreases Ay.  The raw model
    becomes post-peak when this direction reverses.  Retaining the continuous
    segment that contains beta=0 prevents a plotted line from doubling back
    at either end without altering any calculated forces.
    """
    num_delta, num_beta = ay_grid.shape
    center_beta_index = int(np.argmin(np.abs(beta_values)))
    mask = np.zeros_like(ay_grid, dtype=bool)
    tolerance_g = 1e-10

    for i in range(num_delta):
        mask[i, center_beta_index] = True

        # Move left from beta=0: Ay must rise as beta becomes more negative.
        j = center_beta_index
        while j > 0 and ay_grid[i, j - 1] >= ay_grid[i, j] - tolerance_g:
            j -= 1
            mask[i, j] = True

        # Move right from beta=0: Ay must fall as beta becomes more positive.
        j = center_beta_index
        while j < num_beta - 1 and ay_grid[i, j + 1] <= ay_grid[i, j] + tolerance_g:
            j += 1
            mask[i, j] = True

    return mask


def find_trim_points(beta_values, delta_values, ay_grid, mz_grid, valid_mask):
    """Return one near-center Mz=0 steering solution per beta column.

    This deliberately avoids the old behavior of collecting every crossing,
    sorting by Ay, and joining unrelated roots into one artificial line.
    """
    trim_beta = []
    trim_ay = []
    trim_delta = []

    for j, beta_deg in enumerate(beta_values):
        best = None
        for i in range(len(delta_values) - 1):
            if not (valid_mask[i, j] and valid_mask[i + 1, j]):
                continue

            mz_1 = mz_grid[i, j]
            mz_2 = mz_grid[i + 1, j]
            if mz_1 == 0.0:
                fraction = 0.0
            elif mz_1 * mz_2 < 0.0:
                fraction = -mz_1 / (mz_2 - mz_1)
            else:
                continue

            delta_trim = delta_values[i] + fraction * (delta_values[i + 1] - delta_values[i])
            ay_trim = ay_grid[i, j] + fraction * (ay_grid[i + 1, j] - ay_grid[i, j])
            candidate = (abs(delta_trim), beta_deg, ay_trim, delta_trim)
            if best is None or candidate[0] < best[0]:
                best = candidate

        if best is not None:
            _, beta_trim, ay_trim, delta_trim = best
            trim_beta.append(beta_trim)
            trim_ay.append(ay_trim)
            trim_delta.append(delta_trim)

    return np.asarray(trim_beta), np.asarray(trim_ay), np.asarray(trim_delta)


def run_ymd():
    beta_values = np.arange(-SWEEP_LIMIT_DEG, SWEEP_LIMIT_DEG + GRID_STEP_DEG / 2.0, GRID_STEP_DEG)
    delta_values = np.arange(-SWEEP_LIMIT_DEG, SWEEP_LIMIT_DEG + GRID_STEP_DEG / 2.0, GRID_STEP_DEG)

    if DEBUG:
        open("YMD_Debug.txt", "w", encoding="utf-8").close()

    tire = Tire(vp.TireModel, vp.TirePressure_bar)
    shape = (len(delta_values), len(beta_values))
    ay_grid = np.empty(shape)
    mz_grid = np.empty(shape)
    roll_grid = np.empty(shape)
    max_slip_grid = np.empty(shape)
    travel_grid = {
        "FL": np.empty(shape),
        "FR": np.empty(shape),
        "RL": np.empty(shape),
        "RR": np.empty(shape),
    }

    for i, delta_deg in enumerate(delta_values):
        for j, beta_deg in enumerate(beta_values):
            ay_mps2, mz_nm, roll_rad, data = solve(
                VX_MPS,
                beta_deg * DEG2RAD,
                delta_deg * DEG2RAD,
                vp,
                tire,
                DEBUG,
            )
            ay_grid[i, j] = ay_mps2 / vp.Gravity
            mz_grid[i, j] = mz_nm * vp.NM2FTLB
            roll_grid[i, j] = roll_rad * RAD2DEG
            max_slip_grid[i, j] = max(
                abs(data["alpha_FL_deg"]),
                abs(data["alpha_FR_deg"]),
                abs(data["alpha_RL_deg"]),
                abs(data["alpha_RR_deg"]),
            )
            for corner in travel_grid:
                try:
                    travel_grid[corner][i, j] = data[f"Travel_{corner}_mm"]
                except KeyError as error:
                    raise RuntimeError(
                        "This 2DYawMomentDiagram.py requires the matching updated "
                        "YMDSim.py and VehicleParameters.py. Replace all three files "
                        "from YMD_LoopFix together; do not replace only the driver script."
                    ) from error

    branch_mask = monotonic_branch_mask(ay_grid, beta_values)
    if MAX_TIRE_SLIP_DEG is None:
        slip_mask = np.ones_like(branch_mask, dtype=bool)
        limit_text = "monotonic operating branch"
    else:
        slip_mask = max_slip_grid <= MAX_TIRE_SLIP_DEG + 1e-9
        limit_text = f"monotonic branch; tire slip ≤ {MAX_TIRE_SLIP_DEG:g}°"

    valid_mask = branch_mask & slip_mask
    ay_plot = np.where(valid_mask, ay_grid, np.nan)
    mz_plot = np.where(valid_mask, mz_grid, np.nan)
    return beta_values, delta_values, ay_grid, mz_grid, roll_grid, max_slip_grid, travel_grid, valid_mask, ay_plot, mz_plot, limit_text


def plot_ymd(beta_values, delta_values, ay_plot, mz_plot, limit_text):
    plt.figure(figsize=(12, 8))

    for i, delta_deg in enumerate(delta_values):
        if not PLOT_WHOLE_DEGREE_LINES_ONLY or is_whole_degree(delta_deg):
            label = f"δ={int(round(delta_deg))}°" if is_whole_degree(delta_deg) else None
            plt.plot(ay_plot[i, :], mz_plot[i, :], linewidth=1.5, label=label)

    for j, beta_deg in enumerate(beta_values):
        if not PLOT_WHOLE_DEGREE_LINES_ONLY or is_whole_degree(beta_deg):
            plt.plot(ay_plot[:, j], mz_plot[:, j], "--", linewidth=1.2)
            if is_whole_degree(beta_deg) and int(round(beta_deg)) % 2 == 0:
                valid_rows = np.where(np.isfinite(ay_plot[:, j]))[0]
                if len(valid_rows):
                    last = valid_rows[-1]
                    plt.text(ay_plot[last, j], mz_plot[last, j], f"β={int(round(beta_deg))}°", fontsize=8)

    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.axvline(0.0, color="black", linewidth=0.8)
    plt.xlabel("Lateral Acceleration (g)")
    plt.ylabel("Force-Induced Yaw Moment (ft-lbf)")
    plt.title(f"Yaw Moment Diagram — {VX_MPS * vp.MPS2MPH:.2f} mph\n({limit_text})")
    plt.grid(True)
    plt.legend(title="Constant Steering Angle", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()


def plot_trim(trim_beta, trim_ay, trim_delta, limit_text):
    if not len(trim_ay):
        print("No valid Mz=0 trim points inside the requested operating envelope.")
        return

    # Plot in beta order, which preserves one continuous trim branch.
    plt.figure(figsize=(8, 6))
    plt.plot(trim_ay, trim_delta, "o-", markersize=4)
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.axvline(0.0, color="black", linewidth=0.8)
    plt.xlabel("Lateral Acceleration (g)")
    plt.ylabel("Trim Steering Angle δ (deg)")
    plt.title(f"Trim Steering vs Lateral Acceleration\n({limit_text})")
    plt.grid(True)
    plt.tight_layout()


def print_summary(beta_values, delta_values, ay_grid, mz_grid, roll_grid, max_slip_grid, travel_grid, valid_mask, ay_plot, mz_plot):
    i0 = int(np.argmin(np.abs(delta_values)))
    j0 = int(np.argmin(np.abs(beta_values)))
    dbeta_rad = (beta_values[j0 + 1] - beta_values[j0 - 1]) * DEG2RAD
    ddelta_rad = (delta_values[i0 + 1] - delta_values[i0 - 1]) * DEG2RAD
    d_mz_d_beta = -(mz_grid[i0, j0 + 1] - mz_grid[i0, j0 - 1]) / dbeta_rad
    d_mz_d_delta = (mz_grid[i0 + 1, j0] - mz_grid[i0 - 1, j0]) / ddelta_rad

    maximum_index = np.unravel_index(np.nanargmax(ay_plot), ay_plot.shape)
    print(f"Vx = {VX_MPS * vp.MPS2MPH:.2f} mph")
    print("===================================")
    print("YMD CENTER DERIVATIVES")
    print("===================================")
    print(f"dMz/dbeta = {d_mz_d_beta:.2f} ft-lbf/rad")
    print(f"dMz/ddelta = {d_mz_d_delta:.2f} ft-lbf/rad")
    print(f"control-to-stability ratio = {d_mz_d_delta / d_mz_d_beta:.2f}")
    print("===================================")
    print("OPERATING-ENVELOPE LIMITS")
    print("===================================")
    print(f"AY MAX = {np.nanmax(ay_plot):.3f} g")
    print(f"MZ MAX = {np.nanmax(mz_plot):.1f} ft-lbf")
    print(f"Mz at AY max = {mz_plot[maximum_index]:.1f} ft-lbf")
    print(f"valid plotted points = {np.count_nonzero(valid_mask)}/{valid_mask.size}")
    print(f"raw max tire slip = {np.nanmax(max_slip_grid):.2f} deg")
    print(f"raw max roll = {np.nanmax(np.abs(roll_grid)):.3f} deg")
    print("wheel travel at AY max (mm; +compression, -extension):")
    print(
        "  " + ", ".join(
            f"{corner}={travel_grid[corner][maximum_index]:+.2f}"
            for corner in ("FL", "FR", "RL", "RR")
        )
    )

    # Report the true suspension demand across the plotted operating envelope,
    # not just at the single maximum-Ay point.  This excludes the deliberately
    # discarded post-peak branches and any point beyond the tire-slip limit.
    print("maximum valid wheel travel (mm; +compression, -extension):")
    for corner in ("FL", "FR", "RL", "RR"):
        valid_travel = np.where(valid_mask, travel_grid[corner], np.nan)
        compression_index = np.unravel_index(np.nanargmax(valid_travel), valid_travel.shape)
        extension_index = np.unravel_index(np.nanargmin(valid_travel), valid_travel.shape)

        compression_point = (
            f"beta={beta_values[compression_index[1]]:+.2f} deg, "
            f"delta={delta_values[compression_index[0]]:+.2f} deg, "
            f"Ay={ay_grid[compression_index]:+.3f} g"
        )
        extension_point = (
            f"beta={beta_values[extension_index[1]]:+.2f} deg, "
            f"delta={delta_values[extension_index[0]]:+.2f} deg, "
            f"Ay={ay_grid[extension_index]:+.3f} g"
        )
        print(
            f"  {corner}: compression={valid_travel[compression_index]:+.2f} mm "
            f"({compression_point}); extension={valid_travel[extension_index]:+.2f} mm "
            f"({extension_point})"
        )


def main():
    results = run_ymd()
    (
        beta_values,
        delta_values,
        ay_grid,
        mz_grid,
        roll_grid,
        max_slip_grid,
        travel_grid,
        valid_mask,
        ay_plot,
        mz_plot,
        limit_text,
    ) = results

    print_summary(
        beta_values,
        delta_values,
        ay_grid,
        mz_grid,
        roll_grid,
        max_slip_grid,
        travel_grid,
        valid_mask,
        ay_plot,
        mz_plot,
    )
    trim_beta, trim_ay, trim_delta = find_trim_points(
        beta_values, delta_values, ay_grid, mz_grid, valid_mask
    )
    plot_ymd(beta_values, delta_values, ay_plot, mz_plot, limit_text)
    plot_trim(trim_beta, trim_ay, trim_delta, limit_text)
    plt.show()


if __name__ == "__main__":
    main()
