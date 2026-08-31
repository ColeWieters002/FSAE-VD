# FSAE Vehicle Dynamics

Vehicle dynamics simulation for our FSAE car. The goal is to build a model that can predict vehicle behavior from our actual vehicle, suspension, and tire parameters.

## Project Structure

```text
FSAE-VD/
│
├── Python/
│   ├── BicycleSim.py
│   ├── VehicleParameters.py
│   ├── TireFunctions.py
│   └── Tire/
│       └── *.tir
│
└── README.md
```
## File Explination
This section will outline the use and funciton behind each of the models with the inputs and outputs if there are any
### VehicleParameters.py
This file is used to house all of the constants that will be used reguarding the vehicle. Some of the variables that are included stem from the mass of the vehicle, to the trackwidth, to the Yaw Momement of Inertia.  
  
Values yet to be set
* Sprung Roll Inertia
* Sprung Pitch Inertia
* F&R Roll Stiffness
* F&R Heave Stiffness
* Chassis Stiffness
* LLTD
* KPI
* Caster Angle
* F&R Toe

### TireFunctions.py
This file intakes the BicycleSim.py coefficients and translates it into 3 calculators:  
* Fy
  * Converts inputs into Lateral Force
  * Inputs:
    * X = (Slip Angle (deg), Normal Force(lbf), Incination Angle(deg))
    * Lateral Components
    * Nominal Load Coefficients
* Fx
  * Converts inputs into Longitudnal Force
  * Inputs:
    * X = (Slip Angle (deg), Normal Force(lbf), Incination Angle(deg))
    * Aligning Coefficients
    * Lateral Components
    * Nominal Load Coefficients
    * Radius
* Mz
  * Converts inputs into Aligning Force
  * Inputs:
    * X = (Slip Angle (deg), Normal Force(lbf), Incination Angle(deg))
    * Lateral Components
    * Nominal Load Coefficients
### BicycleSim.py
This file creates the tire class and parses the .tir file into sorted arrays for the TireFunctions.py file to use. It also includes the solve function that imports Vehicle Parameters and uses them to solve using scipy's optimize.
* Tire (Class)
  * FY - Outputs lateral force for that tire.
  * FX - Outputs longitudinal force for that tire.
  * Mz - Outputs Aligning force for that tire.
* Solve (Function)
  * Input:
    * Vx - Longitudinal Velocity
    * delta - Steering Angle
    * vp - Vehicle Parameters
    * tire - Tire that is being solved for
  * Output:
    * beta - Slideslip Angle
    * r - Yaw Rate
    * Ay - Lateral Acceleration
### BicycleMain.py
This file combines all previously mentioned files to simulate a steady state estimation for a corner.
* Inputs:
  * Vx - longitudinal Velocity (m/s)
  * delta - Steering Angle (deg)
* Outputs:
  * Sideslip Angle
  * Turn Radius
  * Skidpad Time - Time it takes to do a complete circle

  
