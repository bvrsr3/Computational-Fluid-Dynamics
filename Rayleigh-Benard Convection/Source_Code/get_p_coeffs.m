function [aPp, aE, aW, aN, aS, b] = get_p_coeffs(nx,ny,ust,vst,dxc,dyc,dt)
dx=1/(nx-1);
dy=1/(ny-1);
%dt=0.00001;

  aPp = zeros(nx,ny);
  aE = zeros(nx,ny);  aW = zeros(nx,ny);
  aN = zeros(nx,ny);  aS = zeros(nx,ny);
  b  = zeros(nx,ny);
%   du = zeros(size(ust));  dv = zeros(size(vst));

%   atmp = repmat(dyc', nx,1);
  aE = dy/dx*ones(nx,ny);   %02b = b - ust(2:nx+1,:).*atmp;
  aW = dy/dx*ones(nx,ny);   %b = b + ust(1:nx,:).*atmp;
  %du(2:nx,:) = atmp(2:nx,:)./aPu(2:nx,:);

%   atmp = repmat(dxc, 1, ny);
  aN =dx/dy*ones(nx,ny);   %b = b - vst(:,2:ny+1).*atmp;
  aS = dx/dy*ones(nx,ny);   %b = b + vst(:,1:ny).*atmp;
  %dv(:,2:ny) = atmp(:,2:ny)./aPv(:,2:ny);

  % set boundary coefficients to zero
  aW(1,:) = 0; aE(nx,:) = 0; aS(:,1) = 0; aN(:,ny) = 0;

  aPp = aE + aW + aN + aS;

  for i=1:nx
      for j=1:ny
 
       b(i,j)= -(dy*(ust(i+1,j)-ust(i,j))+dx*(vst(i,j+1)-vst(i,j)))/dt;
          
      end
  end

  % all Neumann boundaries lead to a singular matrix
  % fix this by fixing pressure at (1,1)
  %aPp(1,1) = 1/0.7; 
  aE(1,1) = 0; aW(1,1) = 0; aN(1,1) = 0; aS(1,1) = 0;
  %b(1) = 0;

  %for j=1:nx
  % for i=1:ny
  %   [i j aE(i,j) aW(i,j) aN(i,j) aS(i,j) aPp(i,j)]
  % end
  %end
end
