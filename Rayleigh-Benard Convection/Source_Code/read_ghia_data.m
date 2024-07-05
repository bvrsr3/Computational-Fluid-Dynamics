function [xref, yref, uref, vref] = read_ghia_data(Re)

uall = dlmread('ghiau.dat');
yref = uall(:,1); 
u100 = uall(:,2);
u400 = uall(:,3);
u1k  = uall(:,4);


vall = dlmread('ghiav.dat');
xref = vall(:,1); 
v100 = vall(:,2);
v400 = vall(:,3);
v1k  = vall(:,4);

if(abs(Re-100)<1e-4)
   uref = u100; vref = v100;
elseif(abs(Re-400)<1e-4)
   uref = u400; vref = v400;
elseif(abs(Re-1000)<1e-4)
   uref = u1k; vref = v1k;
else
  'Data not set for Re = ', Re
end

end

%%% ghiau.dat and ghiav.dat contain profiles for the following Re values:
%%% 100, 400, 1000, 3200, 5000, 7500, 10000
