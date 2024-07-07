# Parallelized - 1D Heat diffusion equation 

## Governing equation:
Following is the 1D Heat conduction equation. In this work, the equation is solved numerically using a serial code and parallel code using MPI.
<p align = "center">
$\frac{\partial T}{\partial t} = \kappa \frac{\partial ^2 T}{\partial x^2}$
</p>
where $\kappa$ is the thermal conductivity (W/mK)

A Dirchilet boundary condition is used at the boundaries.
For short times, an analytical solution exists given by $T(x,t) = erf(\frac{x - 0.5}{2\sqrt(t)})$ in this way, the numerical methodology can be validated.

## Serial code:
The governing equation is discretized in space and time. A central Difference scheme ($2^{nd}$ order) is used for space, and the Euler time step ($1^{st}$ order) is used in time. The resulting equation shown below is explicitly updated in time. 
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/7f47be4d-be63-4eb0-a48f-3d0a9bcf4bf3" width = "500" height = "50">
</div>

A figure depicting the node topology is shown below.
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/fb468502-d8e8-43d3-8de5-0bc265c1f28b" width = "600" height = "100">
</div>

## Parallel code and methodology:
The processor topology for the parallel code is shown below:
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/9bc44599-9ee7-4359-b0bc-aa84a4718311" width = "700" height = "125">
</div>
</div>

The following are the main considerations for writing the parallel code:
1. Proc0 receives inputs from the "input file" and broadcasts to other processors
2. Generating global grid and local grid for each processor
3. Enforcing BCs for each processor
4. Display local output
5. Debugging information of each processor
6. Halo exchange of information at processor boundaries
7. Update the equation in each processor
8. Reconstruct all local processor data to proc0 

## Results and Discussion: 
### Validation of the parallel code:
L2 norm: The L2 norm between serial and parallel code at the 100th time step is 0.0 up to the 15th digit place. Shown below is a snippet of the same for a few $x$ values. Hence, the parallel code is now validated.
<div align = "center">
<img src = "https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/41b9ac57-f555-471a-b21d-e50f454d4c22" width = "500" height = "125">
</div>

### Strong and Weak scaling:
Strong scaling is defined as 'Speedup as a function of processor number with fixed problem size'. Weak scaling is defined as 'Speedup as a function of processor number with problem size increasing proportionally to increase in processors'.
                                                  
#### Weak scaling:
                                                      
 ![image](https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/b134f9fb-8dec-49c7-a364-a4af841e1e06)

Comments:
The above is similar to the concept of weak scaling. The number of processors is held constant at 2, and the number of data points to be solved is increased from 10 to 10240. The total time increases linearly as the number of data points for computation increases, and most of the contribution is from the calculation time. The communication time almost remains constant as the number of communication exchanges is fixed.

 
                                                      Strong scaling:
 
![image](https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/395d4190-1752-4774-bdbc-db503193613e)

Comments:
The above result gives an idea of strong scaling. The number of data points to be solved is kept fixed at 10240, and the number of procs is increased from 2 to 8.  Initially, the dominant contribution is from calculation time, but this drops out as the number of processors is increased. The communication time drastically increases beyond 4 procs, and the value keeps oscillating. 
The usual trend is that at larger values of processors, the communication takes up major time as compared to calculation time.
