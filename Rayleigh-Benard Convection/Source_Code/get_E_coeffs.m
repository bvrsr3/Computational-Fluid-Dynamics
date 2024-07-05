function [aTp, aTE, aTW, aTN, aTS, bT] = get_E_coeffs(nx,ny,ust,vst,dxe,dye,dxc,dyc,dt,T)
alpha=1;
dx=1/(nx-1);
dy=1/(ny-1);

  aTp = zeros(nx,ny);
  aTE = zeros(nx,ny);  aTW = zeros(nx,ny);
  aTN = zeros(nx,ny);  aTS = zeros(nx,ny);
  bT  = zeros(nx,ny);

Fw = zeros(size(T)); Fe = zeros(size(T)); Fn = zeros(size(T)); Fs = zeros(size(T));
Dw = zeros(size(T)); De = zeros(size(T)); Dn = zeros(size(T)); Ds = zeros(size(T));
%Sp = zeros(size(u)); Su = zeros(size(u));
bT  = zeros(size(T));

  zo = zeros(size(nx,ny));
  
%cellvol = (dxc(2:nx-1)*dyc(2:ny-1));
% celldy = ones(size(dxe)).*dyc';
% celldx = dxe.*ones(size(dyc'));

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % convective and diffusive fluxes at interiors
Fw(1:nx,2:ny-1) = ust(1:nx,2:ny-1).*dy;          Dw(1:nx,2:ny-1) = repmat(alpha./dxe(1:nx) , 1, ny-1) .*dy;
Fe(1:nx,2:ny-1) = ust(2:nx+1, 2:ny-1).*dy;       De(1:nx,2:ny-1) = repmat(alpha./dxe(1:nx) , 1, ny-2) .*dy;
Fs(1:nx,2:ny-1) = vst(1:nx, 2:ny-1) .*dx;        Ds(1:nx,2:ny-1) = repmat(alpha./dye(1:ny-2),  1,nx).*dx;
Fn(1:nx,2:ny-1) = vst(1:nx,3:ny) .*dx;           Dn(1:nx,2:ny-1) = repmat(alpha./dye(2:ny-1),  1,nx).*dx;

% calculate coefficient values --BCs handled through ghost points
advec_scheme = 1;
if(advec_scheme==1) % upwind scheme
  aTE = De + max(zo, -Fe);  aTN = Dn + max(zo, -Fn);  
  aTW = Dw + max(Fw,  zo);  aTS = Ds + max(Fs,  zo);  
elseif(advec_scheme==2) % C-D scheme
  aTE = De - 0.5*Fe;  aTN = Dn - 0.5*Fn;  
  aTW = Dw + 0.5*Fw;  aTS = Ds + 0.5*Fs;  
elseif(advec_scheme==3) % PowerLaw scheme
  aTE = De .* max(0, (1-0.1*abs(Fe./De)).^5) + max(0, -Fe);
  aTW = Dw .* max(0, (1-0.1*abs(Fw./Dw)).^5) + max(0,  Fw);
  aTN = Dn .* max(0, (1-0.1*abs(Fn./Dn)).^5) + max(0, -Fn);
  aTS = Ds .* max(0, (1-0.1*abs(Fs./Ds)).^5) + max(0, -Fs);
end

 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
bT(1:nx,1:ny) = 0.0; 

  aTp = -(aTE + aTW + aTN + aTS) - (Fe-Fw+Fn-Fs);

  %for j=1:nx
  % for i=1:ny
  %   [i j aE(i,j) aW(i,j) aN(i,j) aS(i,j) aPp(i,j)]
  % end
  %end
end
