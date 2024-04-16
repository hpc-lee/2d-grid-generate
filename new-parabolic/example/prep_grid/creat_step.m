clc;
clear all;
close all;

num_of_step = 21;
dh = -50;
for i=1:num_of_step
 step(i) = dh;
end
% for i=1:10
%   step(i) = -20;
% end
% for i=11:40
%   step(i) = -20-(i-10);
% end
% for i=41:801
%   step(i) = -50;
% end
%creat step file
file_name = '../step_file_2d.txt';
fid=fopen(file_name,'w'); % Output file name 
fprintf(fid,'# number of step\n'); 
fprintf(fid,'%d\n',num_of_step);
fprintf(fid,'# step\n'); 
for i=1:num_of_step
  fprintf(fid,'%.9e \n',step(i));
end
