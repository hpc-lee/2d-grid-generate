% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_km = 1;
flag_topo_z = 1;

nx = 401;

dx = 100;
origin_x = 0;
origin_z = 0;

bz = zeros(nx,2);
for i=1:nx

    bz(i,1) = origin_x + (i-1)*dx;
    bz(i,2) = origin_z;
end

if flag_topo_z
  x0 = 10*1e3;
  x1 = 30*1e3;
  a = 3*1e3;
  H = 6*1e3;
  for i = 1:nx
      x = (i-1)*dx;
      topo1 = H*exp(-(x-x0)^2/a^2) - H*exp(-(x-x1)^2/a^2);
      bz(i,2)=bz(i,2)+topo1;
  end
end

% A=-0.000001;
% [bz]=arc_strech(A,bz);

if flag_printf
    figure(1)   
    if flag_km == 1
      bz1=bz/1e3;
      plot(bz1(:,1),bz1(:,2),'k');
      axis equal;
      xlabel('X axis (km)');
      ylabel('Y axis (km)');
    else 
      bz1=bz;
      plot(bz1(:,1),bz1(:,2),'k');
      axis equal;
      xlabel('X axis (m)');
      ylabel('Y axis (m)');
    end
%     hold on;
%     plot(20,-0.5,'k*');
%     plot(bz1(8/0.05,1),bz1(8/0.05,2),'kv');
%     plot(bz1(12/0.05,1),bz1(12/0.05,2),'kv');
%     plot(bz1(22/0.05,1),bz1(22/0.05,2),'kv');
%     plot(bz1(28/0.05,1),bz1(28/0.05,2),'kv');
%     plot(bz1(32/0.05,1),bz1(32/0.05,2),'kv');
%     text(7,-5,'Vp=3.0 km s{-1}');
%     text(7,-8,'Vs=2.0 km s{-1}');
%     text(7,-11,'rho=1.5 g cm{-3}');
%     set(gcf,'color','white','renderer','painters');
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

