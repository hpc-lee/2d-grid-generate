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

nx = 1001;

p1 = 1;
p2 = 201;
p3 = 301;
p4 = 701;
p5 = 801;
p6 = 1001;

dh = 1;
origin_x = 0;
origin_z = 0;

bz = zeros(nx,2);
% NOTE: x direction tight
% so coord x incremental with x
for i=p1:p2
    bz(i,1) = origin_x + (i-1)*dh;
    bz(i,2) = origin_z;
end
for i=p2+1:p3
    bz(i,1) = origin_x + bz(p2,1) + cos(0.5*pi*90/90)*(i-p2)*dh;
    bz(i,2) = origin_z + bz(p2,2) + sin(0.5*pi*90/90)*(i-p2)*dh;
end
for i=p3+1:p4
    bz(i,1) = origin_x + bz(p3,1) + (i-p3)*dh;
    bz(i,2) = origin_z + bz(p3,2);
end
for i=p4+1:p5
    bz(i,1) = origin_x + bz(p4,1) + cos(0.5*pi*90/90)*(i-p4)*dh;
    bz(i,2) = origin_z + bz(p4,2) - sin(0.5*pi*90/90)*(i-p4)*dh;
end
for i=p5+1:p6
    bz(i,1) = origin_x + bz(p5,1) + (i-p5)*dh;
    bz(i,2) = origin_z + bz(p5,2);
end
% bz(:,2) = smooth(bz(:,2),20);
% bz(:,1) = smooth(bz(:,1),20);
if flag_printf
    figure(1)   
    plot(bz(:,1),bz(:,2),'k');
    hold on;
    plot(200.0, 0.0  , 'kv');
    plot(200.0, 100.0, 'kv');
    plot(300.0, 100.0, 'kv');
    plot(400.0, 100.0, 'kv');
    plot(500.0, 100.0, 'kv');
    plot(600.0, 100.0, 'kv');
    plot(600.0, 0.0  , 'kv');
    text(200.0+10, 0.0  +20,'R1',FontSize=14);
    text(200.0+10, 100.0+20,'R2',FontSize=14);
    text(300.0+10, 100.0+20,'R3',FontSize=14);
    text(400.0+10, 100.0+20,'R4',FontSize=14);
    text(500.0+10, 100.0+20,'R5',FontSize=14);
    text(600.0+10, 100.0+20,'R6',FontSize=14);
    text(600.0+10, 0.0  +20,'R7',FontSize=14);

    plot(350,-100,'k*');
    axis equal;
    set(gcf,'Position',[200,200,650,400]);
    xlabel('X axis (m)',FontSize=15,FontWeight='bold');
    ylabel('Y axis (m)',FontSize=15,FontWeight='bold');
    xlim([-100,900]);
    ylim([-300,200]);
%     title('concave step model',FontSize=15);
    set(gca,'layer','top');
    set(gca,'FontSize',10,FontWeight='bold');
    set(gcf,'color','white','renderer','painters');
    text(0,-50,'Vp=3.0 km s^{-1}',FontSize=14);
    text(0,-100,'Vs=2.0 km s^{-1}',FontSize=14);
    text(0,-150,'\rho=1.5 g cm^{-3}',FontSize=14);
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
   print(gcf,'model4.png','-r300','-dpng');
 end

