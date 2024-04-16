clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which grid profile to plot
subs=[1,1];    
subc=[-1,-1];   % '-1' to plot all points in this dimension
subt=[8,8];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 1;
flag_title  = 0;
scl_daspect = [1 1 1];
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

coordinfo=locate_coord(parfnm,output_dir,subs,subc,subt);
[x,z]=gather_coord(coordinfo,output_dir);

%- set coord unit
if flag_km
   x=x/1e3;
   z=z/1e3;
   str_unit='km';
else
   str_unit='m';
end

%-----------------------------------------------------------
%-- set figure
%-----------------------------------------------------------
hid = figure;
set(hid,'BackingStore','on');

plot(x,z,'k-');
hold on
plot(x',z','k-');
  
xlabel(['X axis (' str_unit ')'],FontSize=15);
ylabel(['Y axis (' str_unit ')'],FontSize=15);
  
set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis equal tight;

% title
if flag_title
    gridtitle='XOY-Grid';
    title(gridtitle);
end

% save and print figure
if flag_print
%     width= 500;
%     height=500;
%     set(gcf,'paperpositionmode','manual');
%     set(gcf,'paperunits','points');
%     set(gcf,'papersize',[width,height]);
%     set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,['gird.png'],'-r300','-dpng');
end
