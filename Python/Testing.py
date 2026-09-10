from BicycleSim import Tire, solve 
import VehicleParameters as vp
import matplotlib.pyplot as plt
import numpy as np


tire = Tire(vp.TireModel, vp.TirePressure_bar)

Vx = 11.75

Fy = []

delta = np.linspace(-13,13,100)

for d in delta:
    Fy.append(solve(Vx, d * vp.DEG2RAD, vp, tire)[2])
    
plt.plot(delta,Fy)
plt.show()
