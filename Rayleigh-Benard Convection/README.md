# Rayleigh-Benard-Convection
## Introduction
Rayleigh-B´enard convection is a type of natural convection flow that is only driven by differences in density due to a temperature gradient. The Rayleigh-B´enard convection occurs in a volume of static fluid in which temperature gradients are introduced. In the case studied, the fluid is kept between two enclosing
parallel plates, and the lower plate is maintained at a higher temperature. The fluid near the lower plate will attain a higher temperature and, therefore, a lower density than the rest of the fluid. Gravity will now force the colder and heavier fluid at the top to sink but is opposed by the viscous forces in the fluid. It is the balance between these two forces that determines if convection will occur or not. If the temperature gradient, and thus the density gradient, is large enough, the gravitational forces will dominate, and instabilities occur. It is the occuring flow patterns that are investigated in this work through numerical methods.

## Governing equations:
The Navier-Stokes, along with the energy equation, governs this problem. 
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/7b4395ad-a990-4996-b145-ce9af7677788" width = "550" height = "350">
</div>

## Boundary conditions:
The geometry and its boundary conditions are shown in the figure below:
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/805d9dac-be8e-455b-8ac6-be2aee78478d" width = "450" height = "350">
</div>

## Numerical Methodology:
A predictor-corrector approach is used to solve the governing equations. The flowchart of the methodology followed is shown below:
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/68be7cc8-ebd9-4d91-b804-1ab887d91dde" width = "450" height = "600">
</div>
