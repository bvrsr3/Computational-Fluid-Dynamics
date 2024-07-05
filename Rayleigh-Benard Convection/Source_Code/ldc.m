clear all
clc
Lx = 1.0; Ly = 1.0;
nx = 30; ny = 30;
dx=1/(nx-1);
dy=1/(ny-1);
% physical parameters
%Re = 400; nu = 1/Re;
%beta=3e-2;
%alpha=25.91e-6;
g=9.81;
Ra=50000;
Pr=10;
%Re=400;
%Re1=100;

Utop=0;
Tbot=1;
K=0.25;

Time_step=2e3;

%numerical parameters
relax_gssor = 0.3;
%relax_p = 0.1; relax_u = 0.7; relax_v = 0.7;
%num_max_simple_iter = 1000;
tolerance = 1e-4;
%tolu = tolerance; tolv = tolerance; 
tolp = tolerance;
%maxiteru = 2; maxiterv = 2;
maxiterp = 1e4;
t=0;

iter=0;
dt=1e-5;

% kdiff=1.798e-5;
% 
% nu=kdiff;
[xe, ye, xc, yc, dxe, dye, dxc, dyc] = set_grid(Lx, Ly, nx, ny, 1);

% set initial conditions
[u, v, p, T] = set_initial_conditions(nx,ny,xe,ye,xc,yc,K);

ust = zeros(size(u)); unew = zeros(size(u));
vst = zeros(size(v)); vnew = zeros(size(v));
pst = zeros(size(p)); pnew = zeros(size(p));
ug=zeros(size(p));vg=zeros(size(p));
p_pr = zeros(size(p));
%T = 0.25*ones(size(T)); 
Tnew = zeros(size(T));

% enforce ghost values according to boundary conditions
[u_gh_we, v_gh_sn, p_gh_we, p_gh_sn] = set_boundary_conditions(nx,ny,u,v,p,T,Utop,Tbot);
u(1,:) = u_gh_we(1,:);   u(nx+1,:) = u_gh_we(2,:);
v(:,1) = v_gh_sn(:,1);   v(:,ny+1) = v_gh_sn(:,2);

T(:,1)=100; T(:,ny)=0.0;

u(:,ny) = 0.0;

% begin iterations
%pst = p; p_pr = p*0;
Tguess = T;

for it = 1:Time_step
    iter=iter+1;

    t=t+dt
 
  %u, v, p

  % enforce ghost values according to boundary conditions
  [u_gh_we, v_gh_sn, p_gh_we, p_gh_sn] = set_boundary_conditions(nx,ny,u,v,p,T,Utop,Tbot);
  u(1,:) = u_gh_we(1,:);   u(nx+1,:) = u_gh_we(2,:);
  v(:,1) = v_gh_sn(:,1);   v(:,ny+1) = v_gh_sn(:,2);

  T(1,:)=T(2,:); T(nx,:)=T(nx-1,:);
  T(:,1)=Tbot;    T(:,ny)=0.0;

  u(:,1) = 0; u(:,ny) = Utop;    % this is necessary to ensure u(1,ny) equals 1.
  v(1,:) = 0; v(nx,:) = 0;
  %u, v, p
  %pst = p;

  % compute u-coefficients
  [aPu, aE, aW, aN, aS, b] = get_u_coeffs(nx,ny,u,v,dxc,dxe,dyc,dye,dt);
  
  for i=2:nx-1
      for j=2:ny-1
        u(i,j)=dt*Pr*(aE(i,j)*u(i+1,j)+aW(i,j)*u(i-1,j)+aN(i,j)*u(i,j+1)+aS(i,j)*u(i,j-1)+(aPu(i,j)+(dx*dy/(dt*Pr)))*u(i,j))/(dx*dy);
      end
 end
   


[u_gh_we, v_gh_sn, p_gh_we, p_gh_sn] = set_boundary_conditions(nx,ny,u,v,p,T,Utop,Tbot);
  u(1,:) = u_gh_we(1,:);   u(nx+1,:) = u_gh_we(2,:);
  v(:,1) = v_gh_sn(:,1);   v(:,ny+1) = v_gh_sn(:,2);
  
  T(1,:)=T(2,:); T(nx,:)=T(nx-1,:);
  T(:,1)=Tbot;    T(:,ny)=0.0;

  u(:,1) = 0; u(:,ny) = Utop;    % this is necessary to ensure u(1,ny) equals 1.
  v(1,:) = 0; v(nx,:) = 0;


  [aPv, aE, aW, aN, aS, b] = get_v_coeffs(nx,ny,u,v,dxc,dxe,dyc,dye,g,Ra,T,Tguess,dt);
  %vst = solve_gssor_v(aPv, aE, aW, aN, aS, b, v, nx, ny, relax_gssor, maxiterv, tolv);
  for i=2:nx-1
      for j=2:ny-1
           v(i,j)=dt*Pr*(aE(i,j)*v(i+1,j)+aW(i,j)*v(i-1,j)+aN(i,j)*v(i,j+1)+aS(i,j)*v(i,j-1)+(aPv(i,j)+(dx*dy/(dt*Pr)))*v(i,j)+b(i,j))/(dx*dy);
      end
  end

  [u_gh_we, v_gh_sn, p_gh_we, p_gh_sn] = set_boundary_conditions(nx,ny,u,v,p,T,Utop,Tbot);
  u(1,:) = u_gh_we(1,:);   u(nx+1,:) = u_gh_we(2,:);
  v(:,1) = v_gh_sn(:,1);   v(:,ny+1) = v_gh_sn(:,2);
  
  T(1,:)=T(2,:); T(nx,:)=T(nx-1,:);
  T(:,1)=Tbot;    T(:,ny)=0.0;

  u(:,1) = 0; u(:,ny) = Utop;    % this is necessary to ensure u(1,ny) equals 1.
  v(1,:) = 0; v(nx,:) = 0;


  % compute p-coefficients
  [aPp, aE, aW, aN, aS, b] = get_p_coeffs(nx,ny,u,v,dxc,dyc,dt);

  p_pr = solve_gssor_p(aPp, aE, aW, aN, aS, b, p_pr, nx, ny, relax_gssor, maxiterp, tolp);
  %p_pr= solve_fullmatrix_p(aPp, aE, aW, aN, aS, b, p_pr, nx, ny);

 [u_gh_we, v_gh_sn, p_gh_we, p_gh_sn] = set_boundary_conditions(nx,ny,u,v,p_pr,T,Utop,Tbot);


  
  % correct u, v and p

  unew(2:nx,2:ny-1) = u(2:nx,2:ny-1)- (p_pr(2:nx,2:ny-1)-p_pr(1:nx-1,2:ny-1))*dt/dx;              %(1-relax_u)*u(2:nx,:) + relax_u*(ust(2:nx,:) + du(2:nx,:).*(p_pr(1:nx-1,:)-p_pr(2:nx,:)));
  vnew(2:nx-1,2:ny) = v(2:nx-1,2:ny)- (p_pr(2:nx-1,2:ny)-p_pr(2:nx-1,1:ny-1))*dt/dy;

  [u_gh_we, v_gh_sn, p_gh_we, p_gh_sn] = set_boundary_conditions(nx,ny,unew,vnew,p_pr,T,Utop,Tbot);

  unew(1,:) = u_gh_we(1,:);   unew(nx+1,:) = u_gh_we(2,:);
  vnew(:,1) = v_gh_sn(:,1);   vnew(:,ny+1) = v_gh_sn(:,2);
  unew(:,1) = 0; unew(:,ny) = Utop;    % this is necessary to ensure u(1,ny) equals 1.
 
  vnew(1,:) = 0; vnew(nx,:) = 0;
  T(1,:)=T(2,:); T(nx,:)=T(nx-1,:);
  T(:,1)=Tbot;   T(:,ny)=0.0;

  
[aTp, aTE, aTW, aTN, aTS, bT] = get_E_coeffs_new(nx,ny,u,v,dxe,dye,dxc,dyc,dt,T);

  for i=2:nx-1
      for j=2:ny-1
  Tnew(i,j)=dt*(aTE(i,j)*T(i+1,j)+aTW(i,j)*T(i-1,j)+aTN(i,j)*T(i,j+1)+aTS(i,j)*T(i,j-1)+(aTp(i,j)+dx*dy/dt)*T(i,j))/(dx*dy);
      end
  end


  [u_gh_we, v_gh_sn, p_gh_we, p_gh_sn] = set_boundary_conditions(nx,ny,unew,vnew,p_pr,Tnew,Utop,Tbot);

  unew(1,:) = u_gh_we(1,:);   unew(nx+1,:) = u_gh_we(2,:);
  vnew(:,1) = v_gh_sn(:,1);   vnew(:,ny+1) = v_gh_sn(:,2);
  unew(:,1) = 0; unew(:,ny) = Utop;    % this is necessary to ensure u(1,ny) equals 1.
  Tnew(:,1) = Tbot; Tnew(:,ny) = 0.0;
  Tnew(1,:)=Tnew(2,:); Tnew(nx,:)=Tnew(nx-1,:);
  
  vnew(1,:) = 0; vnew(nx,:) = 0;

  % prepare for next iteration
  u = unew; v = vnew; p = p_pr; T=Tnew;
VV=zeros(nx,ny);
for i=1:nx
    for j=1:ny

        ug(i,j)=(u(i,j)+u(i+1,j))/2;
        vg(i,j)=(v(i,j)+v(i,j+1))/2;
        VV(i,j)=sqrt(ug(i,j).^2 +vg(i,j).^2);
    end
end

%F=zeros(1e5,1);
if rem(iter,100)==0
% plot contours of T
figure(4), clf
contourf(xe, ye, T(1:nx,1:ny)', 'LineStyle', 'none');
set(gca,'fontsize',14), colorbar, xlabel('x'), ylabel('y'), title(strcat('T-',num2str(nx,'%04d'),+';',num2str(Ra,'%04d'),+';',num2str(Pr,'%04d')))
screen2jpeg(strcat('contours_T_',num2str(nx,'%04d'),'.png'))
%F(iter/100)=im2frame(T);

figure(5), clf
contourf(xe, ye, VV(1:nx,1:ny)', 'LineStyle', 'none');
set(gca,'fontsize',14), colorbar, xlabel('x'), ylabel('y'), title(strcat('V Magnitude-',num2str(nx,'%04d'),+';',num2str(Ra,'%04d'),+';',num2str(Pr,'%04d')))
screen2jpeg(strcat('contours_V_Mag_',num2str(nx,'%04d'),'.png'))
end

end
%VV=zeros(nx,ny);

% for i=1:nx
%     for j=1:ny
% 
%         ug(i,j)=(u(i,j)+u(i+1,j))/2;
%         vg(i,j)=(v(i,j)+v(i,j+1))/2;
%         VV(i,j)=sqrt(ug(i,j).^2 +vg(i,j).^2);
%     end
% end


% plot contours of u
figure(1), clf
contourf(xc, ye, u(2:nx,:)', 'LineStyle', 'none');
set(gca,'fontsize',14), colorbar, xlabel('x'), ylabel('y'), title(strcat('u-',num2str(nx,'%04d'),+';',num2str(Ra,'%04d'),+';',num2str(Pr,'%04d')))
screen2jpeg(strcat('contours_u_',num2str(nx,'%04d'),'.png'))

% plot contours of v
figure(2), clf
contourf(xe, yc, v(:,2:ny)', 'LineStyle', 'none');
set(gca,'fontsize',14), colorbar, xlabel('x'), ylabel('y'), title(strcat('v-',num2str(nx,'%04d'),+';',num2str(Ra,'%04d'),+';',num2str(Pr,'%04d')))
screen2jpeg(strcat('contours_v_',num2str(nx,'%04d'),'.png'))

% plot contours of p
figure(3), clf
contourf(xe, ye, p(1:nx,1:ny)', 'LineStyle', 'none');
set(gca,'fontsize',14), colorbar, xlabel('x'), ylabel('y'), title(strcat('p-',num2str(nx,'%04d'),+';',num2str(Ra,'%04d'),+';',num2str(Pr,'%04d')))
screen2jpeg(strcat('contours_p_',num2str(nx,'%04d'),'.png'))

% figure(5), clf
% contourf(xe, ye, VV(1:nx,1:ny)', 'LineStyle', 'none');
% set(gca,'fontsize',14), colorbar, xlabel('x'), ylabel('y'), title(strcat('V Magnitude-',num2str(nx,'%04d'),+';',num2str(Ra,'%04d'),+';',num2str(Pr,'%04d')))
% screen2jpeg(strcat('contours_V_Mag_',num2str(nx,'%04d'),'.png'))

% % plot contours of T
% figure(4), clf
% contourf(xe, ye, T(1:nx,1:ny)', 'LineStyle', 'none');
% set(gca,'fontsize',14), colorbar, xlabel('x'), ylabel('y'), title(strcat('T-',num2str(nx,'%04d'),+';',num2str(Ra,'%04d'),+';',num2str(Pr,'%04d')))
% screen2jpeg(strcat('contours_T_',num2str(nx,'%04d'),'.png'))

figure(6),clf
quiver(xe,ye,ug,vg);


% [xref, yref, uref, vref] = read_ghia_data(Re);
% [x1ref, y1ref, u1ref, v1ref] = read_ghia_data(Re1);
% 
% 
% %extract u vs y along x=0.5
% [vv, imid] = min(abs(xc-0.5*Lx)); imid = imid+1;
% figure(10), clf
% plot(u(imid,:), ye, 'k-', uref, yref, 'b',u1ref,y1ref,' r');
% set(gca,'fontsize',14), xlabel('u'), ylabel('y'), title(strcat('u-',num2str(nx,'%04d')))
% legend('FV Solver', 'Ghia et al.-400','Ghia et al.-100','Location','SouthEast')
% screen2jpeg(strcat('line_u_y_',num2str(nx,'%04d'),'.png'))
% 
% % extract v vs x along y=0.5
% [vv, jmid] = min(abs(yc-0.5*Ly)); jmid = jmid+1;
% figure(11), clf
% plot(xe, v(:,jmid), 'k-', xref, vref,'b',x1ref,v1ref,'r ');
% set(gca,'fontsize',14), xlabel('x'), ylabel('v'), title(strcat('v-',num2str(nx,'%04d')))
% legend('FV Solver', 'Ghia et al.','Ghia et al.-100','Location','SouthEast')
% screen2jpeg(strcat('line_v_x_',num2str(nx,'%04d'),'.png'))
