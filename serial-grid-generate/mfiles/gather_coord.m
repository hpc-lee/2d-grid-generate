function [x,z] = gather_coord(coordinfo,output_dir)

% check
if ~ exist(output_dir,'dir')
    error([mfilename ': directory ' output_dir ' does not exist']);
end

% load
coordprefix='coord';
nthd=length(coordinfo);
for n=1:nthd
    
    n_i=coordinfo(n).thisid(1); n_k=coordinfo(n).thisid(2);
    i1=coordinfo(n).indxs(1); k1=coordinfo(n).indxs(2);
    i2=coordinfo(n).indxe(1); k2=coordinfo(n).indxe(2);
    subs=coordinfo(n).subs;
    subc=coordinfo(n).subc;
    subt=coordinfo(n).subt;
    fnm_coord=[output_dir,'/',coordprefix,'_px',num2str(n_i),'_pz',num2str(n_k),'.nc'];
    
    if ~ exist(fnm_coord,'file')
       error([mfilename ': file ' fnm_coord 'does not exist']);
    end

    subs=fliplr(subs);subc=fliplr(subc);subt=fliplr(subt);

    %- x coord
    x(k1:k2,i1:i2)=nc_varget(fnm_coord,'x',subs-1,subc,subt);

    %- z coord
    z(k1:k2,i1:i2)=nc_varget(fnm_coord,'z',subs-1,subc,subt);
end

end
