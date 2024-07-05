function [xe, ye, xc, yc, dxe, dye, dxc, dyc] = set_grid(Lx, Ly, nx, ny, type)

  if(type==1) 
    % uniform for now
    dx = Lx/(nx-1); 
    xc = ([1:nx-1]'-0.5)*dx;                               % size is nx-1
    xe = ([0:nx-1]')*dx;                                     % size is nx
    dxc = [2*(xc(1)-xe(1)); diff(xc); 2*(xe(nx)-xc(nx-1))];    % size is nx
    dxe = diff(xe);                                        % size is nx-1

    dy = Ly/(ny-1);
    yc = ([1:ny-1]'-0.5)*dy;                               % size is ny-1
    ye = ([0:ny-1])'*dy;                                     % size is ny
    dyc = [2*(yc(1)-ye(1)); diff(yc); 2*(ye(ny)-yc(ny-1))];    % size is ny
    dye = diff(ye);                                        % size is ny-1
  else
      'Not implemented'
  end


end
