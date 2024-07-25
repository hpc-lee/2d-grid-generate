% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_km = 1;
flag_topo_z = 0;

nx1 = 2501;
num_pml = 00;
nx = nx1 + 2*num_pml; 

dx = 100;
dz = 100;
origin_x = 0;
origin_z = 0;

bz = zeros(nx,2);
for i=1:nx1
    bz(i+num_pml,1) = origin_x + (i-1)*dx;
    bz(i+num_pml,2) = origin_z;
end

if flag_topo_z
  point1 = 5*1e3;
  point2 = 15*1e3;
  L = 4*1e3;
  H = 1.5*1e3;
  for i = 1:nx1
      r1 = sqrt((bz(i+num_pml,1)-point1)^2);
      topo1 = 0;
      if(r1 < L)
          topo1 = H * (1+cos(pi*r1/L));
      end
      bz(i+num_pml,2)=bz(i+num_pml,2)+topo1;

      r2 = sqrt((bz(i+num_pml,1)-point2)^2);
      topo2 = 0;
      if(r2 < L)
          topo2 = H * (1+cos(pi*r2/L));
      end
      bz(i+num_pml,2)=bz(i+num_pml,2)-topo2;
  end
end

[bz] = extend_abs_layer(bz,dx,nx,num_pml);
% A=-0.000001;
% [bz]=arc_strech(A,bz);

if flag_printf
    figure(1)   
    if flag_km == 1
      plot(bz(:,1)/1e3,bz(:,2)/1e3,'k');
      axis equal;
      xlabel('X axis (km)');
      ylabel('Y axis (km)');
    else 
      plot(bz(:,1),bz(:,2),'k');
      axis equal;
      xlabel('X axis (m)');
      ylabel('Y axis (m)');
    end
end

% creat data file
file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# bz coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bz(i,1),bz(i,2));
end
