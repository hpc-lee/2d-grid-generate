% creat 2D boundary data file
% nx nz x1(left) x2(right) z1(bottom) z2(top)

clc;
clear all;
close all;

flag_printf = 1;
flag_topo_z = 1;
flag_km = 1;

nx1 = 801;
nz = 601;
num_pml = 0;
nx = nx1 + 2*num_pml; 

dx = 50;
dz = 50;
origin_x = 0;
origin_z = 0;

bz1 = zeros(nx,2);
bz2 = zeros(nx,2);
for i=1:nx1
    bz1(i+num_pml,1) = origin_x + (i-1)*dx;
    bz1(i+num_pml,2) = origin_z - (nz-1)*dz;

    bz2(i+num_pml,1) = origin_x + (i-1)*dx;
    bz2(i+num_pml,2) = origin_z;
end

if flag_topo_z
  x0 = 10*1e3;
  x1 = 30*1e3;
  a = 3*1e3;
  H = 6*1e3;
  for i = 1:nx
      x = (i-1)*dx;
      topo = H*exp(-(x-x0)^2/a^2) - H*exp(-(x-x1)^2/a^2);
      bz2(i,2)=bz2(i,2)+topo;
  end
end

[bz1] = extend_abs_layer(bz1,dx,nx,num_pml);
[bz2] = extend_abs_layer(bz2,dx,nx,num_pml);

% A=-0.000001;
% [bz1]=arc_strech(A,bz1);
% [bz2]=arc_strech(A,bz2);

if flag_printf
  figure(1)   
  plot(bz1(:,1),bz1(:,2),'k');
  hold on;
  plot(bz2(:,1),bz2(:,2),'k');
  axis equal tight;
  set(gcf,'Position',[200,200,650,400]);
  xlabel('X axis (m)');
  ylabel('Y axis (m)');
  set(gca,'layer','top');
  set(gcf,'color','white','renderer','painters');
end

% creat data file
export_bdry;

% save and print figure
if flag_printf
  print(gcf,'model1.png','-r300','-dpng');
end
