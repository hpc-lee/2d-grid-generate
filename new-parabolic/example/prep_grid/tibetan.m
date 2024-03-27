clc;
clear all;
close all;
topo = load('./topo_srtm.mat');
x = topo.xp';
y = topo.yp';
z = topo.zp';

ny1 = size(x,1);
nx1 = size(x,2);

y_line = 900;
a=2000:2300;
x1 = x(y_line,a);
z1 = z(y_line,a);
z1 = smooth(z1,10);
flag_printf = 1;

num_pml = 0;
n_a = length(a);
nx = n_a + 2*num_pml;

bz = zeros(nx,2);
for i=1:n_a
    bz(i+num_pml,1) = x1(i);
    bz(i+num_pml,2) = z1(i);
end

if flag_printf
    figure(1) 
    plot(bz(:,1),bz(:,2));
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
