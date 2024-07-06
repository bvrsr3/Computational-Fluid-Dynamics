# Problem statement:
The Lid-driven cavity is one of the standard validation test cases for incompressible flow solvers. The Lid-Driven cavity problem has been solved numerically for Reynolds number 100. The code has been written in C language. The experimental data has been obtained from Gaia 1982. The code follows a Finite Volume Framework (FVM). The SIMPLE algorithm is used for the solution. The fluxes at the face are reconstructed using Rhie-Chow interpolation. 

## Governing equations:
The governing equations for this problem are the steady, laminar, incompressible Navier-Stokes equations. 

<p align = "center">
$\frac{\partial}{\partial x}(\rho uu) + \frac{\partial}{\partial y}(\rho vu) = \frac{\partial}{\partial x}\left(\mu \frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left(\mu \frac{\partial u}{\partial y} \right) - \frac{\partial p}{\partial x} + S_u$
</p>

<p align = "center">
$\frac{\partial}{\partial x}(\rho uv) + \frac{\partial}{\partial y}(\rho vv) = \frac{\partial}{\partial x}\left(\mu \frac{\partial v}{\partial x} \right) + \frac{\partial}{\partial y}\left(\mu \frac{\partial v}{\partial y} \right) - \frac{\partial p}{\partial y} + S_v$
</p>

<p align = "center">
$\frac{\partial}{\partial x}(\rho u) + \frac{\partial}{\partial y}(\rho v) = 0$
</p>

## SIMPLE Algorithm:
The SIMPLE algorithm is shown below [1]:
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/4a21c035-a651-4d23-9c33-b5aa7feb6aaa">
</div>

## Results:

<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/11110e79-1fe2-4706-b849-4dd736d62e8b" height = "300" width = "300">
</div>

<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/c3d54888-069d-490b-bc9d-1b723b960c8e" height = "300" width = "300">
</div>

### References:
[1] Versteeg, H. K. "Malalasekera: An Introduction to Computational Fluid Dynamics." Harlow, Essex, England: New York (1995).
