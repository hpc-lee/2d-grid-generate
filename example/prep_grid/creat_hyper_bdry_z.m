clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;

nx1 = 250;
num_pml = 00;
nx = nx1 + 2*num_pml; 

dx = 10;
origin_x = 0;
origin_z = 0;
a_z=0.2;
H_z=0.2*nx*dx;

bdry_z = zeros(nx,2);
% NOTE: x direction tight
% so coord x incremental with x
for i=1:nx1
    bdry_z(i+num_pml,1) = origin_x + (i-1)*dx;
    bdry_z(i+num_pml,2) = origin_z;
end

if  flag_topo_z
    for i = 1:nx1
       bdry_z(i+num_pml,2) = bdry_z(i+num_pml,2) + H_z*exp(-((i-1)/(nx1-1)-0.5)^2/a_z^2);
    end
end

[bdry_z] = extend_abs_layer(bdry_z,dx,nx,num_pml);
A = 0.00001;
[bdry_z] = arc_strech(A,bdry_z);

if flag_printf
    figure(1)   
    plot(bdry_z(:,1),bdry_z(:,2));
    axis equal;
end

% creat data file
file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# bz coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bdry_z(i,1),bdry_z(i,2));
end
