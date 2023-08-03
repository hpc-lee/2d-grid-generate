function [v] = gather_media(parfnm,output_dir,varnm,subs,subc,subt)

% load
fnm_quality=[output_dir,'/',varnm,'.nc'];
    
if ~ exist(fnm_quality,'file')
   error([mfilename ': file ' fnm_quality 'does not exist']);
end

xs = subs(1) -1; 
zs = subs(2) -1; 

xzc = nc_attget(fnm_quality,nc_global,'number_of_points');
xzc = double(xzc);

if(subc(1) == -1)
  xc = floor(xzc(1)/subt(1))-subs(1)+1;
else
  xc = subc(1);
end
if(subc(2) == -1)
  zc = floor(xzc(2)/subt(2))-subs(2)+1;
else
  zc = subc(2);
end
%stride
xt = subt(1);
zt = subt(2);

i1 = 1;
i2 = i1 + xc - 1;
k1 = 1;
k2 = k1 + zc - 1;

v(k1:k2,i1:i2)=nc_varget(fnm_quality,varnm,[zs,xs],[zc,xc],[zt,xt]);

%v=v';

end
