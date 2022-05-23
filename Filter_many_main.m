clear
clc
close all

% filters selected HTC BL data files and stores the batches for training.

%% ======= USER INPUT ==========%%
% Varables to loop over
rus_vec = [2,4,8] ; %ratio of undersampling for each case

% File names
source_file_name_list = {};
input_path_list = {};
source_iter_list = {};
source_case_name_list = {};

    % M5_Re1e4_T02
source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '0,0,0-384,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '385,0,0-768,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '769,0,0-1152,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '1153,0,0-1536,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '1537,0,0-1920,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '1921,0,0-2304,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '2305,0,0-2688,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '2689,0,0-3072,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '3073,0,0-3456,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re1e4_T02';
source_iter_list{end+1} =       '0000260000';
source_file_name_list{end+1} =  '3457,0,0-3841,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

    % M5_Re2e4_T06
source_case_name_list{end+1} =  'M5_Re2e4_T06';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '0,0,0-512,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re2e4_T06';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '513,0,0-1024,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re2e4_T06';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '1025,0,0-1536,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re2e4_T06';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '1537,0,0-2048,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M5_Re2e4_T06';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '2049,0,0-2561,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';


    % M584_Re15e4_T025
source_case_name_list{end+1} =  'M584_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '0,0,0-512,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M584_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '513,0,0-1024,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M584_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '1025,0,0-1536,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M584_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '1537,0,0-2048,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M584_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '2049,0,0-2561,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

    % M6_Re15e4_T025
source_case_name_list{end+1} =  'M6_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '0,0,0-512,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M6_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '513,0,0-1024,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M6_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '1025,0,0-1536,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M6_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '1537,0,0-2048,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M6_Re15e4_T025';
source_iter_list{end+1} =       '0000120000';
source_file_name_list{end+1} =  '2049,0,0-2561,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';


    % M7_Re1e4_T02
source_case_name_list{end+1} =  'M7_Re1e4_T02';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '0,0,0-512,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re1e4_T02';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '513,0,0-1024,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re1e4_T02';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '1025,0,0-1536,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re1e4_T02';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '1537,0,0-2048,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re1e4_T02';
source_iter_list{end+1} =       '0000200000';
source_file_name_list{end+1} =  '2049,0,0-2561,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

    % M7_Re2e4_T03
source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '0,0,0-480,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '481,0,0-960,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '961,0,0-1440,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '1441,0,0-1920,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '1921,0,0-2400,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '2401,0,0-2880,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '2881,0,0-3360,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re2e4_T03';
source_iter_list{end+1} =       '0000100000';
source_file_name_list{end+1} =  '3361,0,0-3841,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';


    % M7_Re15e4_T02
source_case_name_list{end+1} =  'M7_Re15e4_T02';
source_iter_list{end+1} =       '0000180000';
source_file_name_list{end+1} =  '0,0,0-512,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re15e4_T02';
source_iter_list{end+1} =       '0000180000';
source_file_name_list{end+1} =  '513,0,0-1024,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re15e4_T02';
source_iter_list{end+1} =       '0000180000';
source_file_name_list{end+1} =  '1025,0,0-1536,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re15e4_T02';
source_iter_list{end+1} =       '0000180000';
source_file_name_list{end+1} =  '1537,0,0-2048,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M7_Re15e4_T02';
source_iter_list{end+1} =       '0000180000';
source_file_name_list{end+1} =  '2049,0,0-2561,257,447.hdf';
input_path_list{end+1} = '../../raw_data/';

    % M11_Re3e4_T03
source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '0,0,0-480,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '481,0,0-960,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '961,0,0-1440,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '1441,0,0-1920,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '1921,0,0-2400,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '2401,0,0-2880,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '2881,0,0-3360,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';

source_case_name_list{end+1} =  'M11_Re3e4_T03';
source_iter_list{end+1} =       '0000240000';
source_file_name_list{end+1} =  '3361,0,0-3841,257,511.hdf';
input_path_list{end+1} = '../../raw_data/';











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










