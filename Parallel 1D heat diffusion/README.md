                                                       L2 norm:

The L2 norm between serial and parallel code at the 100th time step is 0.0 up to the 15th digit place. Shown below is a snippet of the same.

![image](https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/e643e077-a014-40c5-8a21-d844589341e6)
![image](https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/41b9ac57-f555-471a-b21d-e50f454d4c22)


                                                      Weak scaling:
                                                      
 ![image](https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/b134f9fb-8dec-49c7-a364-a4af841e1e06)

Comments:
The above is similar to the concept of weak scaling. The number of processors is held constant at 2, and the number of data points to be solved is increased from 10 to 10240. The total time increases linearly as the number of data points for computation increases, and most of the contribution is from the calculation time. The communication time almost remains constant as the number of communication exchanges is fixed.

 
                                                      Strong scaling:
 
![image](https://github.com/bvrsr3/Computational-Fluid-Dynamics/assets/137035712/395d4190-1752-4774-bdbc-db503193613e)

Comments:
The above result gives an idea of strong scaling. The number of data points to be solved is kept fixed at 10240, and the number of procs is increased from 2 to 8.  Initially, the dominant contribution is from calculation time, but this drops out as the number of processors is increased. The communication time drastically increases beyond 4 procs, and the value keeps oscillating. 
The usual trend is that at larger values of processors, the communication takes up major time as compared to calculation time.
