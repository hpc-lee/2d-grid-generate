clc;
clear all;
close all;

flag_printf = 1;
topo = importdata('../../../topo_coord.dat');

a=1:500;
x=topo(a,1);
z=topo(a,2);
% smooth topo
% z=smooth(z,20);
figure(1)
plot(x,z,'b');
% hold on;
% plot(x,z2,'r');
% set(gcf,'color','w');
% axis equal;

num_pml = 20;
dx=topo(2,1)-topo(1,1);
dz=dx;

nx1 = length(a);
nx = nx1 + 2*num_pml;

% nz = nx1; % optional
nz = 100;

origin_x = min(x);
origin_z = min(z);
bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
bx1 = zeros(nz,2);
bx2 = zeros(nz,2);
for i=1:nx1
    bz1(i+num_pml,1) = x(i);
    bz1(i+num_pml,2) = origin_z - (nz-1)*dz;

    bz2(i+num_pml,1) = x(i);
    bz2(i+num_pml,2) = z(i);
end

[bz1] = extend_abs_layer(bz1,dx,nx,num_pml);
[bz2] = extend_abs_layer(bz2,dx,nx,num_pml);
dz1 = (bz2(1,2) -bz1(1,2))/(nz-1);
dz2 = (bz2(nx,2)-bz1(nx,2))/(nz-1);

for k=1:nz
    bx1(k,1) = bz1(1,1);
    bx1(k,2) = bz1(1,2) + (k-1) * dz1;

    bx2(k,1) = bz1(nx,1);
    bx2(k,2) = bz1(nx,2) + (k-1) * dz2;
end

if flag_printf
    figure(2)   
    plot(bx1(:,1),bx1(:,2));
    hold on;
    plot(bx2(:,1),bx2(:,2));
    plot(bz1(:,1),bz1(:,2));
    plot(bz2(:,1),bz2(:,2));
    axis equal;
end

% creat data file
export_bdry;
