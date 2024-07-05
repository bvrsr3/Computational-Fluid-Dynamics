function v = solve_gssor_v(aP, aE, aW, aN, aS, b, v, nx, ny, alp_relax, max_iter, tol)

  vnew = zeros(size(v));

  %max_iter = 5; tol = 1e-4;
  for iter=1:max_iter

    % use old values
    vold = v;

    for i=2:nx-1
     for j=2:ny
       v(i,j) = (1-alp_relax)*vold(i,j) + alp_relax*(b(i,j) + aE(i,j)*v(i+1,j) + aW(i,j)*v(i-1,j) + aN(i,j)*v(i,j+1) + aS(i,j)*v(i,j-1))/aP(i,j);
       %if(j==ny-1)
       %[alp_relax i j uold(i,j) b(i,j)   aE(i,j)*u(i+1,j)   aW(i,j)*u(i-1,j)   aN(i,j)*u(i,j+1)   aS(i,j)*u(i,j-1) aP(i,j) u(i,j)]
       %end
     end
    end

    % impose BCs
    v(1,:) = 0.0;          % left wall
    v(nx,:) = 0.0;         % right wall
    v(:,1) = -v(:,2);      % bottom wall
    v(:,ny+1) = -v(:,ny);  % top wall

    % calculate errors and check for convergence
    l2diff = sqrt(sum((v-vold).^2, 'all')/((ny+1)*nx));
    norm_fac = mean(abs(v),'all');
    if(norm_fac < tol) norm_fac = 1; end

    %[iter l2diff/norm_fac]
    %v

    if(l2diff/norm_fac < tol) break; end
  end

  %if(l2diff/norm_fac > tol)  
  %   'Solution not converged. Check details.' 
  %end

end
