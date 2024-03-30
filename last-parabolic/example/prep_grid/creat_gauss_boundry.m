% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;

nx1 = 401;
nz = 201;
num_pml = 00;
nx = nx1 + 2*num_pml;

dx = 100;
dz = 100;
origin_x = 0;
origin_z = 0;

bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
for i=1:nx1
    bz1(i+num_pml,1) = origin_x + (i-1)*dx;
    bz1(i+num_pml,2) = origin_z - (nz-1)*dz;

    bz2(i+num_pml,1) = origin_x + (i-1)*dx;
    bz2(i+num_pml,2) = origin_z;
end

if flag_topo_z
  point1 = 10*1e3;
  point2 = 30*1e3;
  L = 8*1e3;
  H = 3*1e3;
  for i = 1:nx1
      r1 = sqrt((bz2(i+num_pml,1)-point1)^2);
      topo1 = 0;
      if(r1 < L)
          topo1 = H * (1+cos(pi*r1/L));
      end
      bz2(i+num_pml,2)=bz2(i+num_pml,2)+topo1;

      r2 = sqrt((bz2(i+num_pml,1)-point2)^2);
      topo2 = 0;
      if(r2 < L)
          topo2 = H * (1+cos(pi*r2/L));
      end
      bz2(i+num_pml,2)=bz2(i+num_pml,2)-topo2;
  end
end

[bz1] = extend_abs_layer(bz1,dx,nx,num_pml);
[bz2] = extend_abs_layer(bz2,dx,nx,num_pml);

A=-0.000001;
[bz1]=arc_strech(A,bz1);
[bz2]=arc_strech(A,bz2);

if flag_printf
  figure(1)   
  plot(bz1(:,1),bz1(:,2),'k');
  hold on;
  plot(bz2(:,1),bz2(:,2),'k');
  axis equal;
  xlabel('X axis (m)');
  ylabel('Y axis (m)');
end
% creat data file
file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# bz1 coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bz1(i,1),bz1(i,2));
end
fprintf(fid,'# bz1 coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bz2(i,1),bz2(i,2));
end
