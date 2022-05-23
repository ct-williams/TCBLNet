clear
clc
close all

% filters selected HTC BL data files and stores the batches for training.

%% ======= USER INPUT ==========%%
% Varables to loop over
%rus_vec = [2,4,8] ; %ratio of undersampling for each case
rus_vec = [8] ; %ratio of undersampling for each case

% File names
source_file_name_list = {};
input_path_list = {};
source_iter_list = {};
source_case_name_list = {};

source_file_name_list{end+1} = '0,0,0-512,257,447.hdf';
source_iter_list{end+1} = '0000120000';
source_case_name_list{end+1} = 'M584_Re15e4_T025';
input_path_list{end+1} = '../../raw_data/';

source_file_name_list{end+1} = '513,0,0-1024,257,447.hdf';
source_iter_list{end+1} = '0000120000';
source_case_name_list{end+1} = 'M584_Re15e4_T025';
input_path_list{end+1} = '../../raw_data/';

%source_file_name_list{end+1} = '1025,0,0-1536,257,447.hdf';
%source_iter_list{end+1} = '0000120000';
%source_case_name_list{end+1} = 'M584_Re15e4_T025';
%input_path_list{end+1} = '../../raw_data/';
%
%source_file_name_list{end+1} = '1537,0,0-2048,257,447.hdf';
%source_iter_list{end+1} = '0000120000';
%source_case_name_list{end+1} = 'M584_Re15e4_T025';
%input_path_list{end+1} = '../../raw_data/';
%
%source_file_name_list{end+1} = '2049,0,0-2561,257,447.hdf';
%source_iter_list{end+1} = '0000120000';
%source_case_name_list{end+1} = 'M584_Re15e4_T025';
%input_path_list{end+1} = '../../raw_data/';





% Constant variables for all cases
output_filtered_data = 0; % do you want to output filtered data?
filtered_out_path ='../../data_filtered/';

average_flag = 0; % is this an average hdf file?

output_training_data = 1; % do you want to output training data?
training_out_path = '../../training_data/';
mini_batch_size = 256;
xmin = 0.1; % fraction of domain length to ignore at start
xmax = 0.9; % fraction of domain length after which to ignore at end
npoints = 5; % number of wall normal points to train on
randperm_flag = 1; % randomly permute the data to output

% Filter Parameters
%rus =2 ; %ratio of undersampling
number_filters = 5;
number_filters_us = 4;
fil_size = 3.17;





%%=======  Compute  ==========%%

% loop through datafiles
for i = 1:length(source_file_name_list)
    sprintf('Processing file: %s%s/fluid_iter%s/%s \n', input_path_list{i},...
            source_case_name_list{i}, source_iter_list{i},source_file_name_list{i})
    % loop through rus
    for rus = rus_vec
        sprintf('Using rus = %i', rus)
        Filter_HTR_BL_data_func(source_file_name_list{i}, source_iter_list{i},...
                         source_case_name_list{i}, input_path_list{i}, ...
                         output_filtered_data, filtered_out_path, average_flag,...
                         output_training_data, training_out_path, mini_batch_size,...
                         xmin, xmax, npoints, randperm_flag, rus, number_filters,...
                         number_filters_us, fil_size)
    end
end










