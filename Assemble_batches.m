clear
clc
close all

% this function selects the cases you want to train on and then
% generates a bunch of randomized batches from them for training.


data_path_list = {};
%%======= USER INPUT =========%%
data_location = '../../training_data/';

    %---- All of these must be the same for all batches.
mini_batch_size = 256;
num_variables = 242; % (including outputs and features)
    % thes ones are just used for documentation
noutputs = 2; %number of outputs used
npoints = 5; % number of wall normal points used
    %-----

output_path = '../../training_data/selection_1';

% Cases to use
% It is recursive, so it will enter all subdirectories
% and find the batch files to use.
%
% these are the paths on top of the data_location path
data_path_list{end+1} = 'M584_Re15e4_T025';





%%======= Data Processing ===========%%
%find all batches in all data_paths
filelist = [];
for i = 1:length(data_path_list)
    % find all batch files in this path and all sub paths
    rootdir = [data_location,data_path_list{i}];
    filelist_tmp = dir(fullfile(rootdir, '**/*.*'));  %get list of files and folders in any subfolder
    filelist_tmp = filelist_tmp(~[filelist_tmp.isdir]);  %remove folders from list
    filelist_tmp = filelist_tmp(contains({filelist_tmp(:).name}, 'batch_')); % select only batch files

    filelist = [filelist; filelist_tmp];
end


    % create random batch_file, batch row pairs
%    file_row_pairs = zeros(length(filelist)*mini_batch_size,2);
%    ind = 0;
%    for j = 1:length(filelist)
%        for k = 1:mini_batch_size
%            ind = ind + 1;
%            file_row_pairs(ind,:) = [j,k];
%        end
%    end
%    indvec = randperm(length(filelist)*mini_batch_size);
%    file_row_pairs = file_row_pairs(indvec,:);

% read all batches
total_samples = length(filelist)*mini_batch_size;
xz = zeros(total_samples, 2);
indices = zeros(total_samples, 2);
vals = zeros(total_samples,num_variables);
for j = 1:length(filelist)
    start_ind = (j-1)*mini_batch_size+1;
    end_ind = start_ind + mini_batch_size - 1;
    tmpname = [filelist(j).folder, '/', filelist(j).name];
    xz(start_ind:end_ind,:) = h5read(tmpname,'/xz_coords');
    indices(start_ind:end_ind,:) = h5read(tmpname,'/indices');
    vals(start_ind:end_ind,:) = h5read(tmpname,'/data');
end

% reorder into a random order
indvec = randperm(total_samples);
xz = xz(indvec,:);
indices = indices(indvec,:);
vals = vals(indvec,:);

% write to new batch files
if ~exist(output_path, 'dir') 
   mkdir(output_path)
end

for j = 1:length(filelist)
    disp(sprintf('saving batch # %i',j))
    fname_tmp = sprintf('/batch_%0.5i.hdf',j);
    fname = [output_path, fname_tmp];
    % create file
    fileID = H5F.create(fname,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

    % Write info parameters
    h5writeatt(fname,'/','Source file name', ['Combination: ', [data_path_list{:}] ]);
%        h5writeatt(fname,'/','Ratio of undersampling, rus', 0);
%        h5writeatt(fname,'/','number_filters', 0);
%        h5writeatt(fname,'/','number_filters_us', 0);
%        h5writeatt(fname,'/','fil_size', 0);
    h5writeatt(fname,'/','batch number', j);
    h5writeatt(fname,'/','num outputs per set (placed at the start of each row)', noutputs);
    h5writeatt(fname,'/','Number of wall-normal points', npoints);
%        h5writeatt(fname,'/','Fraction of x where data was started', xmin);
%        h5writeatt(fname,'/','Fraction of x where data was ended', xmax);
    h5writeatt(fname,'/','number of features per point used', num_variables - noutputs);
%        h5writeatt(fname,'/','total nx of filtered dataset', 0);
%        h5writeatt(fname,'/','total nz of filtered dataset', 0);
%        h5writeatt(fname,'/','Was the data randomly permuted (0=no)', 1);        
%        h5writeatt(fname,'/','ordering for non-permuted dataset, (x,z): (1,1), (1,2), (1,3)... (2,1), (2,2)...', 1);
%        h5writeatt(fname,'/','dx', dx);
%        h5writeatt(fname,'/','dz', dz);

    % list of variables saved in order
    h5writeatt(fname,'/',['VARIABLE ORDER: ','tau_wall, heatflux_wall, y(1:n), U(1:n), ',...
                          'V(1:n), W(1:n), P(1:n), rho(1:n), T(1:n), ',...
                          'soundspeed(1:n), mu(1:n), kappa(1:n), ',...
                          'Re_y(1:n), (T-Tw)/(Te-Tw) (1:n), ',...
                          'duidxj(1:n)[(1,1), (1,2), (1,3), (2,1),...,(3,3)], ',...
                          'd2uidxjdxk(1:n)[(1,1,1),(1,2,2),(1,3,3),(1,2,1),(1,3,1),',...
                          '(1,2,3),(2,1,1),...,(3,2,3)], ',...
                          'dTdxi(1:n)[(1),(2),(3)], ',...
                          'd2Tdxjdxk(1:n)[(1,1),(2,2),(3,3),(2,1),(3,1),(2,3)].']...
                           , 0);

    % Write some data
        % Initialize the datasets in the hdf file
    h5create(fname,'/data',[mini_batch_size,num_variables]);
    h5create(fname,'/indices',[mini_batch_size,2]);
    h5create(fname,'/xz_coords',[mini_batch_size,2]);

        % write the batch
    start_ind = (j-1)*mini_batch_size+1;
    end_ind = start_ind + mini_batch_size - 1;

    h5write(fname, '/data',vals(start_ind:end_ind,:),[1,1],[mini_batch_size,num_variables]);
    h5write(fname, '/indices',indices(start_ind:end_ind,:),[1,1],[mini_batch_size,2]);
    h5write(fname, '/xz_coords',xz(start_ind:end_ind,:),[1,1],[mini_batch_size,2]);

    % Close the hdf5 file
    H5F.close(fileID) 
end



return

% testing
figure
hold on
vals = [];
xz = [];
indices = [];
for i = 8:11
    tmpname = sprintf('./batch_%0.5i.hdf',i);
    xz = [xz; h5read(tmpname,'/xz_coords')];
    indices = [indices; h5read(tmpname,'/indices')];
    vals = [vals; h5read(tmpname,'/data')];
end
x = linspace(min(xz(:,1)),max(xz(:,1)),500);
z = linspace(min(xz(:,2)),max(xz(:,2)),200);
[xq,yq] = meshgrid(x, z); % interpolation points grid
vq = griddata(xz(:,1),xz(:,2),vals(:,9),xq,yq,'natural'); %U
%     vq = griddata(xz(:,1),xz(:,2),vals(:,2),xq,yq,'natural');
Z = squeeze(vq);
contourf(xq,yq,Z);

figure;
plot(xz(:,1),xz(:,2),'x')

figure;
plot(indices(:,1),indices(:,2),'x')




