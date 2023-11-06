clc;
clear all;
close all;

% p1        p2        p5            p6
% ...........         ...............
%           .         .
%           .         .
%           ...........    
%           p3        p4

flag_printf = 1;

nx = 250;

p1 = 1;
p2 = 70;
p3 = 100;
p4 = 150;
p5 = 180;
p6 = 250;

dh = 10;
origin_x = 0;
origin_z = 0;

bz = zeros(nx,2);
% NOTE: x direction tight
% so coord x incremental with x
for i=p1:p2
    bz(i,1) = origin_x + (i-1)*dh;
    bz(i,2) = origin_z;
end
for i=p2+1:p3
    bz(i,1) = origin_x + bz(p2,1) + cos(pi*100/200)*(i-p2)*dh;
    bz(i,2) = origin_z + bz(p2,2) - sin(pi*100/200)*(i-p2)*dh;
end
for i=p3+1:p4
    bz(i,1) = origin_x + bz(p3,1) + (i-p3)*dh;
    bz(i,2) = origin_z + bz(p3,2);
end
for i=p4+1:p5
    bz(i,1) = origin_x + bz(p4,1) + cos(pi*100/200)*(i-p4)*dh;
    bz(i,2) = origin_z + bz(p4,2) + sin(pi*100/200)*(i-p4)*dh;
end
for i=p5+1:p6
    bz(i,1) = origin_x + bz(p5,1) + (i-p5)*dh;
    bz(i,2) = origin_z + bz(p5,2);
end
bz(:,2) = smooth(bz(:,2),10);
bz(:,1) = smooth(bz(:,1),10);
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
