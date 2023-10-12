function [x,z] = gather_coord(output_dir,subs,subc,subt)

% load
fnm_coord=[output_dir,'/','coord_px0_pz0.nc'];
    
if ~ exist(fnm_coord,'file')
   error([mfilename ': file ' fnm_coord 'does not exist']);
end

xzc = nc_attget(fnm_coord,nc_global,'count_of_physical_points');
xzc = double(xzc);

xt = subt(1);
if(subc(1) == -1)
  xc = ceil((xzc(1)-subs(1)+1)/xt);
else
  xc = subc(1);
end
xs = subs(1)-1;

zt = subt(2);
if(subc(2) == -1)
  zc = ceil((xzc(2)-subs(2)+1)/zt);
else
  zc = subc(2);
end
zs = subs(2)-1;

i1 = 1;
i2 = i1 + xc - 1;
k1 = 1;
k2 = k1 + zc - 1;

% check dimension size of x,y,z to detmine cart or curv
xvar_info = ncinfo(fnm_coord,'x');
num_dim_x = length(xvar_info.Dimensions);

x(k1:k2,i1:i2)=nc_varget(fnm_coord,'x',[zs,xs],[zc,xc],[zt,xt]);

%- z coord
zvar_info = ncinfo(fnm_coord,'z');
num_dim_z = length(zvar_info.Dimensions);
z(k1:k2,i1:i2)=nc_varget(fnm_coord,'z',[zs,xs],[zc,xc],[zt,xt]);
