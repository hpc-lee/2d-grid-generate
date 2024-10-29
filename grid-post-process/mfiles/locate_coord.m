function coordinfo = locate_coord(parfnm,output_dir,subs,subc,subt)

% check parameter file exist
if ~ exist(parfnm,'file')
    error([mfilename ': file ' parfnm ' does not exist']);
end

% read parameters file
jsontext=fileread(parfnm);
par=jsondecode(jsontext);
num_of_grid = length(par.input_grids_info);

total_nx = 0;
total_nz = 0;
if(strcmp(par.merge_direction,"x") == 1)
  for i=1:num_of_grid
    total_nx = total_nx + par.input_grids_info(i).number_of_grid_points(1);
    total_nz = par.input_grids_info(1).number_of_grid_points(2);
  end
  total_nx = total_nx - num_of_grid + 1;
  ngik=[total_nx,total_nz];
end


if(strcmp(par.merge_direction,"z") == 1)
  for i=1:num_of_grid
    total_nz = total_nz + par.input_grids_info(i).number_of_grid_points(2);
    total_nx = par.input_grids_info(1).number_of_grid_points(1);
  end
  total_nz = total_nz - num_of_grid + 1;
  ngik=[total_nx,total_nz];
end

if(par.flag_sample == 1)
  ngik(1) = ngik(1)*par.sample_factor_xi;
  ngik(2) = ngik(2)*par.sample_factor_zt;
end

gsubs = subs;
gsubt = subt;
gsubc = subc;

% reset count=-1 to total number
indx=find(subc==-1);
gsubc(indx)=ceil((ngik(indx)-gsubs(indx)+1)./gsubt(indx));

gsube=gsubs+(gsubc-1).*gsubt;

% search the nc file headers to locate the threads/processors
coordprefix='coord';
coordlist=dir([output_dir,'/',coordprefix,'*.nc']);
n=0;
for i=1:length(coordlist)
    
    coordnm=[output_dir,'/',coordlist(i).name];
    xzs=double(ncreadatt(coordnm, '/', 'global_index_of_first_physical_points'));
    xs=xzs(1);
    zs=xzs(2);
    xzc=double(ncreadatt(coordnm, '/', 'count_of_physical_points'));
    xc=xzc(1);
    zc=xzc(2);
    xarray=[xs:xs+xc-1];
    zarray=[zs:zs+zc-1];
    if (length(find(xarray>=gsubs(1)-1 & xarray<=gsube(1)-1)) ~= 0 && ...
        length(find(zarray>=gsubs(2)-1 & zarray<=gsube(2)-1)) ~= 0)
        n=n+1;
        px(n)=str2num(coordlist(i).name( strfind(coordlist(i).name,'px' )+2 : ...
                                         strfind(coordlist(i).name,'_pz')-1));
        pz(n)=str2num(coordlist(i).name( strfind(coordlist(i).name,'pz' )+2 : ...
                                         strfind(coordlist(i).name,'.nc')-1));
    end
    
end

nthd=0;
for ip=1:length(px)
    
    nthd=nthd+1;
    
    coordnm=[output_dir,'/',coordprefix,'_px',num2str(px(ip)),...
            '_pz',num2str(pz(ip)),'.nc'];
    xzs=double(ncreadatt(coordnm, '/', 'global_index_of_first_physical_points'));
    xs=xzs(1);
    zs=xzs(2);
    xzc=double(ncreadatt(coordnm, '/', 'count_of_physical_points'));
    xc=xzc(1);
    zc=xzc(2);
    xe=xs+xc-1;
    ze=zs+zc-1;
    
    gxarray=gsubs(1):gsubt(1):gsube(1);
    gxarray=gxarray-1;
    gzarray=gsubs(2):gsubt(2):gsube(2);
    gzarray=gzarray-1;
    
    coordinfo(nthd).thisid=[px(ip),pz(ip)];
    i=find(gxarray>=xs & gxarray<=xe);
    k=find(gzarray>=zs & gzarray<=ze);
    coordinfo(nthd).indxs=[i(1),k(1)];
    coordinfo(nthd).indxe=[i(end),k(end)];
    coordinfo(nthd).indxc=coordinfo(nthd).indxe-coordinfo(nthd).indxs+1;
    
    coordinfo(nthd).subs=[ gxarray(i(1))-xs+1, ...
                           gzarray(k(1))-zs+1 ];
    coordinfo(nthd).subc=coordinfo(nthd).indxc;
    coordinfo(nthd).subt=gsubt;
    
    coordinfo(nthd).fnmprefix=coordprefix;
    
end


end

