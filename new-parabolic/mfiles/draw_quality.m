clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project4/test.json';
output_dir='../project4/output';

% which grid profile to plot
% subs=[150,701];     % % index 1:nx 1:nz
% subc=[100,100];   % '-1' to plot all points in this dimension
% subt=[1,1];
subs=[1,1];     % % index 1:nx 1:nz
subc=[200,-1];   % '-1' to plot all points in this dimension
subt=[2,1];
% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 0;
scl_daspect = [1 1 1];
clrmp       = 'parula';

% varable to plot 
%  'orth', 'jacobi', 'ratio', 'smooth_xi', 
%  'smooth_zt', 'step_xi', 'step_zt'
% varnm = 'jacobi';
varnm = 'orth';
% varnm = 'smooth_zt';
%-----------------------------------------------------------
%-- load coord
%-----------------------------------------------------------

qualityinfo=locate_quality(parfnm,output_dir,subs,subc,subt);
[x,z]=gather_coord(qualityinfo,output_dir);

v=gather_quality(qualityinfo,output_dir,varnm);

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

% figure plot
hid=figure;
set(hid,'BackingStore','on');
x1 = flipdim(x,1);
z1 = flipdim(z,1);
v1 = flipdim(v,1);
pcolor(x1,z1,v1);

xlabel(['X axis (' str_unit ')']);
ylabel(['Y axis (' str_unit ')']);

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% shading
% shading interp;
% shading flat;
% colorbar range/scale

min_v = min(min(v));
max_v = max(max(v));
% caxis([min_v,max_v]);

% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis equal;
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
end
% title
if flag_title
    title(varnm,'interpreter','none');
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[varnm '.png'],'-dpng');
end
