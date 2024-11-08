function [v] = gather_quality(qualityinfo,output_dir,varnm)

% check
if ~ exist(output_dir,'dir')
    error([mfilename ': directory ' output_dir ' does not exist']);
end

% load
varprefix=varnm;
nthd=length(qualityinfo);
for n=1:nthd
    
    n_i=qualityinfo(n).thisid(1); n_k=qualityinfo(n).thisid(2);
    i1=qualityinfo(n).indxs(1); k1=qualityinfo(n).indxs(2);
    i2=qualityinfo(n).indxe(1); k2=qualityinfo(n).indxe(2);
    subs=qualityinfo(n).subs;
    subc=qualityinfo(n).subc;
    subt=qualityinfo(n).subt;
    fnm_var=[output_dir,'/',varprefix,'_px',num2str(n_i),'_pz',num2str(n_k),'.nc'];
    
    if ~ exist(fnm_var,'file')
       error([mfilename ': file ' fnm_var 'does not exist']);
    end

    %- vaule 
    v(i1:i2,k1:k2)=ncread(fnm_var,varnm,subs,subc,subt);

end

end
