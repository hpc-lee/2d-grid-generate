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

nx1 = length(a);
nx = nx1 + 2*num_pml;

origin_x = min(x);
origin_z = min(z);
bz = zeros(nx,2);
for i=1:nx1
    bz(i+num_pml,1) = x(i);
    bz(i+num_pml,2) = z(i);
end

[bz] = extend_abs_layer(bz,dx,nx,num_pml);

if flag_printf
    figure(2)   
    plot(bz(:,1),bz(:,2));
%     axis equal;
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
