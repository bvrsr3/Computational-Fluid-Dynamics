function p = solve_gssor_p(aP, aE, aW, aN, aS, b, p, nx, ny, alp_relax, max_iter, tol)

  ppad = zeros(nx+2, ny+2);
  ppad(2:nx+1,2:ny+1) = p;

  %aP, aE, aW, aN, aS, b, p, ppad

  %max_iter = 10000; tol = 1e-4;
  for iter=1:max_iter

    % set ghost cells in ppad for homog. Neumann conditions
    ppad(1,:) = ppad(2,:);    ppad(nx+2,:) = ppad(nx+1,:);
    ppad(:,1) = ppad(:,2);    ppad(:,ny+2) = ppad(:,ny+1);

    % use old values
    ppadold = ppad;

    for i=1:nx
     for j=1:ny
       ipad = i+1; jpad=j+1;
       ppad(ipad,jpad) = (1-alp_relax)*ppadold(ipad,jpad) + alp_relax*(b(i,j) + aE(i,j)*ppad(ipad+1,jpad) + aW(i,j)*ppad(ipad-1,jpad) + aN(i,j)*ppad(ipad,jpad+1) + aS(i,j)*ppad(ipad,jpad-1))/aP(i,j);
       %[i j uold(i,j) b(i,j)   aE(i,j)*u(i+1,j)   aW(i,j)*u(i-1,j)   aN(i,j)*u(i,j+1)   aS(i,j)*u(i,j-1) aP(i,j) u(i,j)]
     end
    end

    % calculate errors and check for convergence
    l2diff = sqrt(sum((ppad-ppadold).^2, 'all')/((nx+2)*(ny+2)));
    norm_fac = mean(abs(ppad),'all');
    if(norm_fac < tol) 
        norm_fac = 1; 
    
   end

    %[iter l2diff/norm_fac]
    %u
    %ppadold, ppad

    if(l2diff/norm_fac < tol) 
        break; 
    end
  end

  if(l2diff/norm_fac > tol)  
     'Pressure Solution not converged. Check details.'
%      strcat('iter = ', num2str(iter), '; l2diff = ', num2str(l2diff), '; norm_fac = ', num2str(norm_fac))
%      strcat('norm_l2diff = ', num2str(l2diff/norm_fac), '; tol = ', num2str(tol))
  end
  p = ppad(2:nx+1,2:ny+1);

end
