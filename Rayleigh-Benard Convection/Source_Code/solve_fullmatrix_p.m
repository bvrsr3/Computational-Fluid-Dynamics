function un2  = solve_fullmatrix_p(aP, aE, aW, aN, aS, b, u, nx, ny)

  aP, aE, aW, aN, aS, b, u%, ppad
  ntot = nx*ny;

  afull = zeros(ntot); bfull = zeros(ntot,1);

  % interior points
  for i=2:nx-1
    for j=2:ny-1
      indP = (j-1)*nx + i;
      indN = (j-0)*nx + i;        indS = (j-2)*nx + i;
      indE = (j-1)*nx + i+1;      indW = (j-1)*nx + i-1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j); 
      afull(indP, indE) = -aE(i,j);      afull(indP, indW) = -aW(i,j);
      afull(indP, indN) = -aN(i,j);      afull(indP, indS) = -aS(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % west boundary
  for i=1
    for j=2:ny-1
      indP = (j-1)*nx + i;
      indN = (j-0)*nx + i;        indS = (j-2)*nx + i;
      indE = (j-1)*nx + i+1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j); 
      afull(indP, indE) = -aE(i,j);
      afull(indP, indN) = -aN(i,j);      afull(indP, indS) = -aS(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % east boundary
  for i=nx
    for j=2:ny-1
      indP = (j-1)*nx + i;
      indN = (j-0)*nx + i;        indS = (j-2)*nx + i;
                                  indW = (j-1)*nx + i-1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j);  
                                        afull(indP, indW) = -aW(i,j);
      afull(indP, indN) = -aN(i,j);      afull(indP, indS) = -aS(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % south boundary
  for i=2:nx-1
    for j=1
      indP = (j-1)*nx + i;
      indN = (j-0)*nx + i;
      indE = (j-1)*nx + i+1;      indW = (j-1)*nx + i-1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j); %+ aS(i,j)*u(i,j-1);
      afull(indP, indE) = -aE(i,j);     afull(indP, indW) = -aW(i,j);
      afull(indP, indN) = -aN(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % north boundary
  for i=2:nx-1
    for j=ny
      indP = (j-1)*nx + i;
                                  indS = (j-2)*nx + i;
      indE = (j-1)*nx + i+1;      indW = (j-1)*nx + i-1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j); %+aN(i,j)*u(i,j+1);  
      afull(indP, indE) = -aE(i,j);      afull(indP, indW) = -aW(i,j);
                                        afull(indP, indS) = -aS(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % NW boundary
  for i=1
    for j=ny
      indP = (j-1)*nx + i;
                                  indS = (j-2)*nx + i;
      indE = (j-1)*nx + i+1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j);  
      afull(indP, indE) = -aE(i,j);
                                        afull(indP, indS) = -aS(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % NE boundary
  for i=nx
    for j=ny
      indP = (j-1)*nx + i;
                                  indS = (j-2)*nx + i;
                                  indW = (j-1)*nx + i-1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j);  
                                        afull(indP, indW) = -aW(i,j);
                                        afull(indP, indS) = -aS(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % SE boundary
  for i=nx
    for j=1
      indP = (j-1)*nx + i;
      indN = (j-0)*nx + i;  
                                  indW = (j-1)*nx + i-1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j);  
                                        afull(indP, indW) = -aW(i,j);
      afull(indP, indN) = -aN(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  % SW boundary
  for i=1
    for j=1
      indP = (j-1)*nx + i;
      indN = (j-0)*nx + i; 
      indE = (j-1)*nx + i+1;

      afull(indP, indP) = aP(i,j);      bfull(indP) = b(i,j);  
      afull(indP, indE) = -aE(i,j);
      afull(indP, indN) = -aN(i,j);

      % safeguard
      indP = 0; indE = 0; indW = 0; indN = 0; indS = 0;
    end
  end

  %max(afull,[],'all')
  %afull = 0.7*afull;
  afull, bfull
  %cond(afull)
  %writematrix(afull,'afull.dat')

  ufull = afull\bfull;

  % fill un2 matrix
  un2 = u;
  for i=1:nx
    for j=1:ny
      indP = (j-1)*nx + i;
      un2(i,j) = ufull(indP);
    end
  end
  %un2
  %u
  %bfull - afull*ufull
end

