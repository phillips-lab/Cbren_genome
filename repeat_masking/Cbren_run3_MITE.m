% Test example for detectMITE.
clear;clc;

data_file = '/projects/phillipslab/ateterina/Cbren/final_repeats/Cbren.genome.fasta';
genome_name = 'CBREN4';

tic;
	do_MITE_detection(data_file,'-genome',genome_name,'-cpu',24)
runtime = toc;

fid = fopen('detectMITE.Runtime.txt','a');
fprintf(fid,'-------%s-------\n',genome_name);
fprintf(fid,'Runtime: %f s\n',runtime);
fclose(fid);
exit;
