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

nx = 801;

p1 = 1;
p2 = 201;
p3 = 301;
p4 = 501;
p5 = 601;
p6 = 801;

dh = 1;
origin_x = 0;
origin_z = 0;

bz = zeros(nx,3);
% NOTE: x direction tight
% so coord x incremental with x
for i=p1:p2
    bz(i,1) = origin_x + (i-1)*dh;
    bz(i,2) = origin_z;
    bz(i,3) = 0;
end
for i=p2+1:p3
    bz(i,1) = origin_x + bz(p2,1) + cos(0.5*pi*90/90)*(i-p2)*dh;
    bz(i,2) = origin_z + bz(p2,2) - sin(0.5*pi*90/90)*(i-p2)*dh;
    bz(i,3) = 0;
end
for i=p3+1:p4
    bz(i,1) = origin_x + bz(p3,1) + (i-p3)*dh;
    bz(i,2) = origin_z + bz(p3,2);
    bz(i,3) = 0;
end
for i=p4+1:p5
    bz(i,1) = origin_x + bz(p4,1) + cos(0.5*pi*90/90)*(i-p4)*dh;
    bz(i,2) = origin_z + bz(p4,2) + sin(0.5*pi*90/90)*(i-p4)*dh;
    bz(i,3) = 0;
end
for i=p5+1:p6
    bz(i,1) = origin_x + bz(p5,1) + (i-p5)*dh;
    bz(i,2) = origin_z + bz(p5,2);
    bz(i,3) = 0;
end

if flag_printf
    figure(1)   
    plot(bz(:,1),bz(:,2)+500,'k');
    axis equal;
end

filename='model3';
igesout(bz,filename);
