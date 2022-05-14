clear
% close all
clc

% Routine to plot P, rho, T, u,v,w from an HDF file.

%% File name
% filename = './data_tar/0,0,0-512,257,447.hdf';
filename = './data_tar/fluid_iter0000120000/513,0,0-1024,257,447.hdf';
% filename = './data_tar/ZAverages/0,0,0-512,257,0.hdf';
% filename = './data_filtered/0,0,0-512,257,447_filtered.hdf';
% filename = './data_filtered/0,0,0-512,257,447_filtered_rus8_numfil5_numfilus4_filsz3.17.hdf';

%% Variables to plot
% lists the available variables
% % tmp = h5info(filename);
% % tmp.Datasets.Name 

% Instantaneous options:
varnames = {'pressure', 'rho', 'temperature', 'velocity'};
% varnames = {'velocity'};

% Some of the average options:
% varnames = {'mu_avg', 'Ma', 'velocity_avg','temperature_avg','pressure_avg'};
% disp('WARNING, the x and y locations in the average file are normalized
% somehow')

% Some of the filtered options
% varnames = {'pressure', 'rho', 'temperature', 'velocity_x', 'velocity_y', 'velocity_z', 'mu', 'soundspeed', 'du1dx2', 'd2u(2)dxidxj(5)'};


file_type = 1; % 1=DNS data, 2=Average data, 3=filtered data


%% User inputs
% xind = 250; % xlocation for the yz plane
% yind = 30;  % ylocation for the xz plane
% zind = 224; % zlocation for the xy plane

% plane_type = 'xy'; % 'xy', 'xz', 'yz'
% index = 28;%28;%112; %224;

plane_type = 'xz'; % 'xy', 'xz', 'yz'
index = 30;%4; %30;

% plane_type = 'yz'; % 'xy', 'xz', 'yz'
% index = 250;

%% Read and plot
% Get grid size from centerCoordinates
if file_type == 3
    tmp = h5info(filename,'/centerCoordinates_x');
else
    tmp = h5info(filename,'/centerCoordinates');
end
datasize = tmp.Dataspace.Size;
nx = datasize(1);
ny = datasize(2);
nz = datasize(3);

if file_type == 3
    x = h5read(filename,'/centerCoordinates_x');
    y = h5read(filename,'/centerCoordinates_y');
    z = h5read(filename,'/centerCoordinates_z');
    xyz(1,:,:,:) = x;
    xyz(2,:,:,:) = y;
    xyz(3,:,:,:) = z;
    clear x y z
else
    xyz = h5read(filename,'/centerCoordinates');
    if file_type ==2
        xyz(1,:,:) = repmat(xyz(1,:,end)',1,ny);
    end
end

for i = 1:length(varnames)
    disp(varnames{i})
    
    tmp = h5info(filename,['/',varnames{i}]);
    idim = tmp.Datatype.Size/8; % Assumes all data are double floats
    
    % Read plane of the data to plot
    if strcmp(plane_type,'xy')
        % xy plane
        data_plane = h5read(filename,['/',varnames{i}],[1,1,index], [nx,ny,1],[1,1,1]);
        xplot = squeeze(xyz(1,:,:,index))';
        yplot = squeeze(xyz(2,:,:,index))';
        xlabelstr = 'x';
        ylabelstr = 'y';
    elseif strcmp(plane_type,'xz')
        % xz plane
        data_plane = h5read(filename,['/',varnames{i}],[1,index,1], [nx,1,nz],[1,1,1]);
        xplot = squeeze(xyz(1,:,index,:))';
        yplot = squeeze(xyz(3,:,index,:))';
        xlabelstr = 'x';
        ylabelstr = 'z';
    elseif strcmp(plane_type,'yz')
        % yz plane
        data_plane = h5read(filename,['/',varnames{i}],[index,1,1], [1,ny,nz],[1,1,1]);
        xplot = squeeze(xyz(3,index,:,:))';
        yplot = squeeze(xyz(2,index,:,:))';
        xlabelstr = 'z';
        ylabelstr = 'y';
    end
    
    % plot them
    for l = 1:idim
        figure
%         if length(size(data_plane))==3
        if idim>1
            varplot = squeeze(data_plane(l,:,:))';
        else
            varplot = squeeze(data_plane)';
        end
        pcolor(xplot,yplot,varplot)
%         contourf(xplot,yplot,varplot)
        shading flat
        colorbar
        title(sprintf('%s, %i, at index %i',varnames{i}, l, index))
        xlabel(xlabelstr)
        ylabel(ylabelstr)
        axis equal
        axis tight
        adfa = 1;
        
    end

end




