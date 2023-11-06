clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;

nx1 = 300;
num_pml = 00;
nx = nx1 + 2*num_pml; 

dx = 10;
origin_x = 0;
origin_z = 0;
a_z=0.2;
H_z=0.3*nx*dx;

bz = zeros(nx,2);
% NOTE: x direction tight
% so coord x incremental with x
for i=1:nx1
    bz(i+num_pml,1) = origin_x + (i-1)*dx;
    bz(i+num_pml,2) = origin_z;
end

if  flag_topo_z
    for i = 1:nx1
       bz(i+num_pml,2) = bz(i+num_pml,2) + H_z*exp(-((i-1)/(nx1-1)-0.5)^2/a_z^2);
    end
end

[bz] = extend_abs_layer(bz,dx,nx,num_pml);
A = 0.00001;
[bz] = arc_strech(A,bz);

if flag_printf
    figure(1)   
    plot(bz(:,1),bz(:,2));
    axis equal;
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
