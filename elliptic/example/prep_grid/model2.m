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


p1 = 1;
p2 = 201;
p3 = 301;
p4 = 501;
p5 = 601;
p6 = 801;

dh = 1;
origin_x = 0;
origin_z = 0;

nx = 801;
nz = 401;
bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
bx1 = zeros(nz,2);
bx2 = zeros(nz,2);
% NOTE: x direction tight
% so coord x incremental with x
for i=p1:p2
    bz2(i,1) = origin_x + (i-1)*dh;
    bz2(i,2) = origin_z;
end
for i=p2+1:p3
    bz2(i,1) = origin_x + bz2(p2,1) + cos(0.5*pi*77/90)*(i-p2)*dh;
    bz2(i,2) = origin_z + bz2(p2,2) - sin(0.5*pi*77/90)*(i-p2)*dh;
end
for i=p3+1:p4
    bz2(i,1) = origin_x + bz2(p3,1) + (i-p3)*dh;
    bz2(i,2) = origin_z + bz2(p3,2);
end
for i=p4+1:p6
    bz2(i,1) = origin_x + bz2(p4,1) + cos(0.5*pi*77/90)*(i-p4)*dh;
    bz2(i,2) = origin_z + bz2(p4,2) + sin(0.5*pi*77/90)*(i-p4)*dh;
end
for i=p5+1:p6
    bz2(i,1) = origin_x + bz2(p5,1) + (i-p5)*dh;
    bz2(i,2) = origin_z + bz2(p5,2);
end

for i=1:801
  bz1(i,1) = bz2(i,1);
  bz1(i,2) = origin_z - (nz-1)*dh;
end


dz1 = (bz2(1,2)-bz1(1,2))/(nz-1);
dz2 = (bz2(nx,2)-bz1(nx,2))/(nz-1);

for k=1:nz
    bx1(k,1) = bz1(1,1);
    bx1(k,2) = bz1(1,2) + (k-1)*dz1;

    bx2(k,1) = bz1(nx,1);
    bx2(k,2) = bz1(nx,2) + (k-1)*dz2;
end


% bz(:,2) = smooth(bz(:,2),20);
% bz(:,1) = smooth(bz(:,1),20);

if flag_printf
    figure(1)   
    plot(bx1(:,1),bx1(:,2));
    hold on;
    plot(bx2(:,1),bx2(:,2));
    plot(bz1(:,1),bz1(:,2));
    plot(bz2(:,1),bz2(:,2));
    axis equal tight;
end
% creat data file
export_bdry;
