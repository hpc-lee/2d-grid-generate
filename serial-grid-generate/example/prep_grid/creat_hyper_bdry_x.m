clc;
clear all;
close all;

flag_printf = 1;
flag_topo_x = 1;

nz = 250;

dz = 10;
origin_x = 0;
origin_z = 0;
a_x=0.2;
H_x=0.2*nz*dz;

bx = zeros(nz,2);
% NOTE: z direction upward
% so coord z incremental with k
for k=1:nz
    bx(k,1) = origin_x;
    bx(k,2) = origin_z + (k-1)*dz;
end

if  flag_topo_x
    for k = 1:nz
       bx(k,1) = bx(k,1) + H_x*exp(-((k-1)/(nz-1)-0.5)^2/a_x^2);
    end
end

A = 0.00001;
[bx] = arc_strech(A,bx);

if flag_printf
    figure(1)   
    plot(bx(:,1),bx(:,2));
    axis equal;
end

% creat data file
file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nz number\n'); 
fprintf(fid,'%d\n',nz);
fprintf(fid,'# bx coords\n'); 
for i=1:nz
  fprintf(fid,'%.9e %.9e\n',bx(i,1),bx(i,2));
end
