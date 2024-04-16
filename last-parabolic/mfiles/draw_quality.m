clear all;
close all;
clc;
addmypath

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which grid profile to plot
subs=[1,1];     % % index 1:nx 1:nz
subc=[-1,-1];   % '-1' to plot all points in this dimension
subt=[2,2];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 1;
flag_clb    = 1;
flag_title  = 1;
scl_daspect = [1 1 1];
clrmp       = 'parula';

% varable to plot 
%  'orth', 'jacobi', 'ratio', 'smooth_xi', 
%  'smooth_zt', 'step_xi', 'step_zt'
% varnm = 'jacobi';
% varnm = 'orth';
varnm = 'smooth_xi';
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

pcolor(x,z,v);

xlabel(['X axis (' str_unit ')'],FontSize=15);
ylabel(['Y axis (' str_unit ')'],FontSize=15);

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');
set(gcf,'Position',[200,200,650,400]);
% shading
shading interp;
% shading flat;
% colorbar range/scale
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis equal tight;
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
end
% title
if flag_title
%     title('Orthogonality',FontSize=15);
      title('Smooth\_xi',FontSize=15);
end

% save and print figure
if flag_print
%   width= 500;
%   height=500;
%   set(gcf,'paperpositionmode','manual');
%   set(gcf,'paperunits','points');
%   set(gcf,'papersize',[width,height]);
%   set(gcf,'paperposition',[0,0,width,height]);
  print(gcf,[varnm,'.png'],'-r300','-dpng');
end
