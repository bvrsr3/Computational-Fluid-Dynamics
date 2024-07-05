function u = solve_gssor_u(aP, aE, aW, aN, aS, b, u, nx, ny, alp_relax, max_iter, tol)

  unew = zeros(size(u));

  %max_iter = 5; tol = 1e-4;
  for iter=1:max_iter

    % use old values
    uold = u;

    %u

    for i=2:nx
     for j=2:ny-1
       u(i,j) = (1-alp_relax)*uold(i,j) + alp_relax*(b(i,j) + aE(i,j)*u(i+1,j) + aW(i,j)*u(i-1,j) + aN(i,j)*u(i,j+1) + aS(i,j)*u(i,j-1))/aP(i,j);
       %if(j==ny-1)
       %[alp_relax i j uold(i,j) b(i,j)   aE(i,j)*u(i+1,j)   aW(i,j)*u(i-1,j)   aN(i,j)*u(i,j+1)   aS(i,j)*u(i,j-1) aP(i,j) u(i,j)]
       %end
     end
    end

    % impose BCs
    u(1,:) = -u(2,:); % left wall
    u(nx+1,:) = -u(nx,:); % right wall
    u(1:nx+1,1) = 0.0;      % bottom wall
    u(1:nx+1,ny) = 1.0;     % top

    % calculate errors and check for convergence
    l2diff = sqrt(sum((u-uold).^2, 'all')/((nx+1)*ny));
    norm_fac = mean(abs(u),'all');
    if(norm_fac < tol) norm_fac = 1; end

    %[iter l2diff/norm_fac]
    %u

    if(l2diff/norm_fac < tol) break; end
  end

  %if(l2diff/norm_fac > tol)  
  %   'Solution not converged. Check details.' 
  %end

end
