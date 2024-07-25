clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project1/test.json';
output_dir='../project1/output';

% which grid profile to plot
subs=[1,1];    
subc=[-1,-1];   % '-1' to plot all points in this dimension
subt=[8,8];

% figure control parameters
flag_km     = 0;
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
% 
plot(x,z,'k-');
hold on
plot(x',z','k-');
  
xlabel(['X axis (' str_unit ')']);
ylabel(['Z axis (' str_unit ')']);
  

xlabel(['X axis (' str_unit ')'],FontSize=15);
ylabel(['Y axis (' str_unit ')'],FontSize=15);

set(gca,'layer','top');
set(gca,'FontSize',10,FontWeight='bold');
set(gcf,'color','white','renderer','painters');
set(gcf,'Position',[200,200,650,400]);
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis equal tight;
% text(-100,-20,'a)',FontSize=20);

% title
if flag_title
    gridtitle='XOZ-Grid';
    title(gridtitle);
end

% save and print figure
if flag_print
    print(gcf,['grid.png'],'-dpng');
end
