% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_x = 1;
flag_topo_z = 1;

nx1 = 300;
nz = 300;
num_pml = 00;
nx = nx1 + 2*num_pml; 

dx = 10;
dz = 10;
origin_x = 0;
origin_z = 0;

bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
bx1 = zeros(nz,2);
bx2 = zeros(nz,2);
for i=1:nx1
    bz1(i+num_pml,1) = origin_x + (i-1)*dx;
    bz1(i+num_pml,2) = origin_z - (nz-1)*dz;

    bz2(i+num_pml,1) = origin_x + (i-1)*dx;
    bz2(i+num_pml,2) = origin_z;
end

if flag_topo_z
  point = origin_x + floor(nx1/2)*dx;
  L = 0.4*(nx-1)*dx;
  H = 0.2*(nz-1)*dz;
  for i = 1:nx1
      r1 = sqrt((bz2(i+num_pml,1)-point)^2);
      topo = 0;
      if(r1 < L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      bz2(i+num_pml,2)=bz2(i+num_pml,2)+topo;
      bz1(i+num_pml,2)=bz1(i+num_pml,2)-topo;
  end
end

[bz1] = extend_abs_layer(bz1,dx,nx,num_pml);
[bz2] = extend_abs_layer(bz2,dx,nx,num_pml);

dz1 = (bz2(1,2)-bz1(1,2))/(nz-1);
dz2 = (bz2(nx,2)-bz1(nx,2))/(nz-1);

for k=1:nz
    bx1(k,1) = bz1(1,1);
    bx1(k,2) = bz1(1,2) + (k-1)*dz1;

    bx2(k,1) = bz1(nx,1);
    bx2(k,2) = bz1(nx,2) + (k-1)*dz2;
end

if flag_topo_x
  point = origin_z - floor(nz/2)*dz;
  L = 0.4*(nz-1)*dz;
  H = 0.2*(nx-1)*dx;
  for k = 1:nz
      r1 = sqrt((bx2(k,2)-point)^2);
      topo = 0;
      if(r1 < L)
          topo = 0.5*H * (1+cos(pi*r1/L));
      end
      bx2(k,1)=bx2(k,1)+topo;
      bx1(k,1)=bx1(k,1)-topo;
  end
end

A=-1;
[bx1]=arc_strech(A,bx1);
[bx2]=arc_strech(A,bx2);

if flag_printf
    figure(1)   
    plot(bx1(:,1),bx1(:,2));
    hold on;
    plot(bx2(:,1),bx2(:,2));
    plot(bz1(:,1),bz1(:,2));
    plot(bz2(:,1),bz2(:,2));
    axis equal;
end
% creat data file
export_bdry;
