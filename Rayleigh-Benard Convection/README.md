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

## Results:
The working fluid is air (Prandtl number, $Pr$ = 0.71). The Rayleigh number $(Ra)$ is increased from 10 to 10,000. At $Ra$ = 10, the temperature difference is not sufficient to overcome the viscous forces and induce motion in the fluid. At $Ra$ = 10,000, a convection roll gets established.

The following figures describe the contours for a 50 X 50 grid, $Ra$ = 10, $Pr$ = 0.71
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/09d47fa1-7f62-4381-b7f8-4627812fd8f6" width = "450" height = "200">  <img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/354890e9-e265-42ae-8ad4-f6e60d022cfc" width = "450" height = "200">
</div>

The following figures describe the contours for a 50 X 50 grid, $Ra$ = 10,000, $Pr$ = 0.71
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/11f15254-fa56-4087-a162-6c6d8eddf54e" width = "450" height = "200">  <img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/50604856-1e68-48f5-a2c3-76362f841d09" width = "450" height = "200">
</div>

## Discussion:
1. It is observed that there exists a critical Rayleigh number below which conduction heat transfer is dominant over convective heat transfer.
2. Beyond the critical Rayleigh number convective mode of heat transfer is dominant. A clockwise convection current is due to numerical directionality.
3. The bands of the temperature contour represent Isotherms. For increased Rayleigh number, these isotherms indicate a strong convective nature of the flow, which further promotes mixing, as observed.

