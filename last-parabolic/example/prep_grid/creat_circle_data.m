% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

nx = 200;
r1 = 50;
r2 = 250;
printf = 1;

for i=1:nx
    bz1(i,1) = r1*cos(-(i-1)/(nx-1)*2*pi);
    bz1(i,2) = r1*sin(-(i-1)/(nx-1)*2*pi);

    bz2(i,1) = r2*cos(-(i-1)/(nx-1)*2*pi);
    bz2(i,2) = 0.5*r2*sin(-(i-1)/(nx-1)*2*pi);
end

if printf == 1
    figure(1)   
    plot(bz1(:,1),bz1(:,2));
    hold on;
    plot(bz2(:,1),bz2(:,2));
    axis equal;
end

% creat data file
export_bdry;
