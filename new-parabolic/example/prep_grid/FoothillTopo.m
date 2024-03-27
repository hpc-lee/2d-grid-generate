close all; 
clear all;
clc; 

addmypath;

flag_printf=1;
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

% a=800:1100; % important topo zone
% a=900:1100;
a=1:1668;
x = x(a);
z = z(a);
z = smooth(z,10);

% plot(x(a),z(a),'r');
% hold on;
% plot(x(a),z1(a),'b');
% axis equal;

nx = length(a);

origin_x = min(x);
origin_z = min(z);

bz = zeros(nx,2);
for i=1:nx
    bz(i,1) = x(i);
    bz(i,2) = z(i);
end

A = 0.00001;
[bz] = arc_strech(A,bz);

if flag_printf
    figure(1) 
    plot(bz(:,1)/1e3,bz(:,2)/1e3);
    axis equal;
    set(gcf,'color','w');
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
