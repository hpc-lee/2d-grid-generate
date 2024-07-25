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
p2 = 301;
p3 = 601;

dh = 1;
origin_x = 0;
origin_z = 0;

nx = 601;
nz = 401;
bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
bx1 = zeros(nz,2);
bx2 = zeros(nz,2);
% NOTE: x direction tight
% so coord x incremental with x
for i=p1:p2
    bz2(i,1) = origin_x + (i-p1)*dh;
    bz2(i,2) = origin_z;
end
for i=p2+1:p3
    bz2(i,1) = bz2(p2,1) + cos(0.5*pi*90/90)*(i-p2)*dh;
    bz2(i,2) = bz2(p2,2) + sin(0.5*pi*90/90)*(i-p2)*dh;
end

for i=1:nx
  bz1(i,1) = bz2(i,1) + (nz-1)*dh;
  bz1(i,2) = bz2(i,2) - (nz-1)*dh;
end


dx1 = (bz2(1,1)-bz1(1,1))/(nz-1);
dx2 = (bz2(nx,1)-bz1(nx,1))/(nz-1);
dz1 = (bz2(1,2)-bz1(1,2))/(nz-1);
dz2 = (bz2(nx,2)-bz1(nx,2))/(nz-1);

for k=1:nz
    bx1(k,1) = bz1(1,1) + (k-1)*dx1;
    bx1(k,2) = bz1(1,2) + (k-1)*dz1;

    bx2(k,1) = bz1(nx,1) + (k-1)*dx2;
    bx2(k,2) = bz1(nx,2) + (k-1)*dz2;
end


% bz(:,2) = smooth(bz(:,2),20);
% bz(:,1) = smooth(bz(:,1),20);

if flag_printf
    figure(1)   
    plot(bx1(:,1),bx1(:,2),'r');
    hold on;
    plot(bx2(:,1),bx2(:,2),'r--');
    plot(bz1(:,1),bz1(:,2),'b');
    plot(bz2(:,1),bz2(:,2),'b--');
    axis equal tight;
end
% creat data file
export_bdry;
