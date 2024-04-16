file_name = '../data_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# nx number\n'); 
fprintf(fid,'%d\n',nx);
fprintf(fid,'# bz1 coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bz1(i,1),bz1(i,2));
end
fprintf(fid,'# bz1 coords\n'); 
for i=1:nx
  fprintf(fid,'%.9e %.9e\n',bz2(i,1),bz2(i,2));
end
