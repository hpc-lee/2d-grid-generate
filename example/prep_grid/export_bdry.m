% creat data file
file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%g\n',nx);
fprintf(fid,'# nz number\n'); 
fprintf(fid,'%g\n',nz);

fprintf(fid,'# bx1 coords\n'); 
for k=1:nz
  fprintf(fid,'%g %g\n',bx1(k,1),bx1(k,2));
end

fprintf(fid,'# bx2 coords\n'); 
for k=1:nz
  fprintf(fid,'%g %g\n',bx2(k,1),bx2(k,2));
end

fprintf(fid,'# bz1 coords\n'); 
for i=1:nx
  fprintf(fid,'%g %g\n',bz1(i,1),bz1(i,2));
end

fprintf(fid,'# bz2 coords\n'); 
for i=1:nx
  fprintf(fid,'%g %g\n',bz2(i,1),bz2(i,2));
end
