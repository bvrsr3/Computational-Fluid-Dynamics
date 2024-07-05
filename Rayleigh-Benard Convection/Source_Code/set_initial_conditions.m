function [u, v, p, T] = set_initial_conditions(nx,ny,xe,ye,xc,yc,K)

  u = zeros(nx+1,ny);
  v = zeros(nx,ny+1);
  p = zeros(nx,ny);
%T = zeros(nx,ny);
T=K*ones(nx,ny);
end
