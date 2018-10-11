% Sequence Alignment Program
% Author: Alaric Wu:yxw1242
% Assignment 2, problem 5

clear all
delimiter = ' ';
inputfile = fopen('input.dat');
C = textscan(inputfile,'%s\n%d %d %d\n%s\n%s');
fclose(inputfile);

[score,number_of_op_solu,op_align] = do_alignment(C{1}{1},C{5}{1},C{6}{1},C{2},C{3},C{4});
global table_m