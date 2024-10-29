clear all;
close all;
clc;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which grid profile to plot
subs=[1,1];     % % index 1:nx 1:nz
subc=[-1,-1];   % '-1' to plot all points in this dimension
subt=[1,1];

% figure control parameters
flag_km     = 0;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
clrmp       = 'parula';
% scl_caxis = [20 90];
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
set(gca,'FontSize',10,FontWeight='bold');
set(gcf,'color','white');
set(gcf,'Position',[200,200,650,400]);
shading interp;
% coorbar range/scale
if exist('scl_caxis')
    caxis(scl_caxis);
end
axis equal;
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
%     set(get(cid,'Title'),'string','degree');
end
if flag_title
%     title('Orthogonality',FontSize=15);
    title('Smooth\_\xi',FontSize=15);
%     title('Smooth\_\eta',FontSize=15);
end

% save and print figure
if flag_print
  print(gcf,[varnm,'.png'],'-r400','-dpng');
end


