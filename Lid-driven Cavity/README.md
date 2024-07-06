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
### Validation
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/ba20fc54-4f89-4b5d-a0ee-a0a805ccc468" height = "300" width = "300">    <img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/67bc0f97-e486-4b24-a4a2-261ddcb1c676" height = "300" width = "300">
</div>

### Contours
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/a4fab3c0-c163-45c7-b0b4-67595a1c960c" height = "300" width = "300">  <img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/b9cc586b-c7a0-4b69-8766-e6545fbcbb8b" height = "300" width = "300">  <img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/3a58092e-750d-4a9a-a941-cbeb07ca32e4">
</div>



### References:
[1] Versteeg, H. K. "Malalasekera: An Introduction to Computational Fluid Dynamics." Harlow, Essex, England: New York (1995).
