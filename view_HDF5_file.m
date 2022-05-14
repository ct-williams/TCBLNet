clear
close all
clc

% Routine to inspect and plot an HDF file to see what is in it.

% filename = '/home/timflint/Documents/Extra_projects/CME230_Machine_learning/data_tar/0,0,0-512,257,447.hdf';
filename = './data_tar/0,0,0-512,257,447.hdf';
% filename = './data_tar/ZAverages/0,0,0-512,257,0.hdf';


% print info
infostruct = h5info(filename)

% print all data info
h5disp(filename)

% Read them
datainfo = infostruct.Datasets;
datanames = {datainfo.Name};
datatypes = {datainfo.Datatype};
dataspaces = {datainfo.Dataspace};

for i = 1:length(datanames)
    disp(datanames{i})
    
    % get grid size
    datasize = dataspaces{i}.Size;
    nx = datasize(1);
    ny = datasize(2);
    nz = datasize(3);
    
    % get individual data size
    ni = datatypes{i}.Size/8; % Assumes all data are double floats
    
    % Read xy plane of the data to plot
    datai = h5read(filename,['/',datanames{i}],[1,1,1], [nx,ny,1],[1,1,1]);
    
    for l = 1:ni
        figure
        if length(size(datai))==3
            pcolor(squeeze(datai(l,:,:))')
        else
            pcolor(datai')
        end
        shading flat
        colorbar
        title(sprintf('%s, %i',datanames{i}, l))
        xlabel('x index')
        ylabel('y index')
    end

end

%i=4 is a probelm


