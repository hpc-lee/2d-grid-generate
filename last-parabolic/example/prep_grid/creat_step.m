clc;
clear all;
close all;

nz = 401;
num_of_step = nz-1;
for i=1:num_of_step
  step(i) = 1;
end
% for i=1:10
%  step(i) = -20;
% end
% for i=11:40
%  step(i) = -20-(i-10);
% end
% for i=41:num_of_step
%  step(i) = -50;
% end
sum_step = sum(step);

% normalization step 
step_nor=step/sum_step;
sum_step_nor = sum(step_nor);
if((sum_step_nor-1)>1e-8)
  error("step set is error, pelase check and reset");
end
% creat step file
file_name = '../step_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# number of step\n'); 
fprintf(fid,'%d\n',num_of_step);
fprintf(fid,'# step\n'); 
for i=1:num_of_step
  fprintf(fid,'%.9e \n',step_nor(i));
end
