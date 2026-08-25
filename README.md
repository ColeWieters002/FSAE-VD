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
This file deciphers the .tir file that stores all of the tire data and translates it into 3 calculators:  
* Fy
  * Inputs:
    * Slip Angle
    * 

* Fx
* Mz
With these values 
