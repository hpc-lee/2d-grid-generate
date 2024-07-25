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

nx = 601;

p1 = 1;
p2 = 301;
p3 = 601;

dh = 1;
origin_x = 0;
origin_z = 0;

bz = zeros(nx,2);
% NOTE: x direction tight
% so coord x incremental with x
for i=p1:p2
    bz(i,1) = origin_x + (i-p1)*dh;
    bz(i,2) = origin_z;
end
for i=p2+1:p3
    bz(i,1) = bz(p2,1) + cos(0.5*pi*90/90)*(i-p2)*dh;  
    bz(i,2) = bz(p2,2) + sin(0.5*pi*90/90)*(i-p2)*dh;
end
% bz(:,2) = smooth(bz(:,2),20);
% bz(:,1) = smooth(bz(:,1),20);
if flag_printf
    figure(1)   
    plot(bz(:,1),bz(:,2),'k');
    axis equal tight;
    xlabel('X axis (m)');
    ylabel('Y axis (m)');
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

if flag_printf
   print(gcf,'model2.png','-r300','-dpng');
 end

