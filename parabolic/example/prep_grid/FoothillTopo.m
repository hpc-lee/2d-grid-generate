close all; 
clear all;
clc; 


flag_printf=1;
num_pml = 20;
%
[Vp,SegyTraceH,SegyH]=ReadSegy('../../../../velocity.segy');
[row, col] = size(Vp);

dz=-10;
dx=15;
% pcolor(Vp(1:5:row,1:5:col));

% obtain topography
for i = 1 : col
    for j = 1:row
        a = Vp(j,i);
        if (abs(a-4000.0)>0.1)
            break;
        end
    end
    z(i)=j*dz;
    x(i)=i*dx;
end

z = smooth(z,20);

nx1 = length(z);
nx = nx1 + 2*num_pml;

nz = 601; % optional

% origin_x = min(x);
% origin_z = min(z);

bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
for i=1:nx1
    bz1(i+num_pml,1) = x(i);
    bz1(i+num_pml,2) = (nz-1)*dz;

    bz2(i+num_pml,1) = x(i);
    bz2(i+num_pml,2) = z(i);
end

[bz1] = extend_abs_layer(bz1,dx,nx,num_pml);
[bz2] = extend_abs_layer(bz2,dx,nx,num_pml);

if flag_printf
    figure(1) 
    plot(bz1(:,1),bz1(:,2));
    hold on;
    plot(bz2(:,1),bz2(:,2));
    axis equal;
    set(gcf,'color','w');
end

% creat data file
export_bdry;




