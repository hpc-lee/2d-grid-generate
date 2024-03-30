clc;
clear all;
close all;

nz = 201; 
num_of_step = nz-1;
for i=1:num_of_step
  step(i) = 1/num_of_step;
end
sum_step = sum(step);
if((sum_step-1)>1e-8)
  error("step set is error, pelase check and reset");
end
% for i=1:10
%   step(i) = -1;
% end
% for i=11:50
%   step(i) = -5;
% end
% for i=51:100
%   step(i) = -10;
% end
% for i=101:200
%   step(i) = -20;
% end
% creat step file
file_name = '../step_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# number of step\n'); 
fprintf(fid,'%d\n',num_of_step);
fprintf(fid,'# step\n'); 
for i=1:num_of_step
  fprintf(fid,'%.9e \n',step(i));
end
