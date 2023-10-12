function [v] = gather_quality(output_dir,varnm,subs,subc,subt)

% load
fnm_quality=[output_dir,'/',varnm,'_px0_pz0.nc'];
    
if ~ exist(fnm_quality,'file')
   error([mfilename ': file ' fnm_quality 'does not exist']);
end

xzc = nc_attget(fnm_quality,nc_global,'count_of_physical_points');
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

v(k1:k2,i1:i2)=nc_varget(fnm_quality,varnm,[zs,xs],[zc,xc],[zt,xt]);

%v=v';

end