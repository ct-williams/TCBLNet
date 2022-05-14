close all;
clear;
clc
% NOTE: HTR advances the dimensional equations.


%% File names
source_file_name = '513,0,0-1024,257,447.hdf';
% source_file_name = '513,0,0-1024,257,447.hdf';
input_path = '../data_tar/fluid_iter0000120000/';

filtered_out_path ='../data_filtered/';

average_flag = 0; % is this an average hdf file?

output_traning_data = 1; % do you want to output training data?
training_out_path = '../training_data/';
mini_batch_size = 256;
xmin = 0.1; % fraction of domain length to ignore at start
xmax = 0.9; % fraction of domain length after which to ignore at end
npoints = 5; % number of wall normal points to train on
randperm_flag = 1; % randomly permute the data to output

%% Filter Parameters
rus =4 ; %ratio of undersampling
number_filters = 5;
number_filters_us = 4;
fil_size = 3.17;


%% Setup

% set file names
fname_in = [ input_path, source_file_name ];
tmpstr = sprintf('_filtered_rus%i_numfil%i_numfilus%i_filsz%0.2f',rus,number_filters,number_filters_us,fil_size);
filtered_file_name = [ source_file_name(1:end-4), tmpstr,'.hdf'];
fname_out = [ filtered_out_path, filtered_file_name ];

% Get grid size from centerCoordinates
tmp = h5info(fname_in,'/centerCoordinates');
datasize = tmp.Dataspace.Size;
nxg = datasize(1);
nyg = datasize(2);
nzg = datasize(3);

% undersampled size
% nx = nxg/rus;
nx = floor((nxg-1)/rus);
% ny = 1+(nyg-1)/rus;
ny = floor(nyg/rus);
nz = floor(nzg/rus);

% Get coordinates and undersample at same time
xyz = h5read(fname_in,'/centerCoordinates',[2,1,1], [nx,ny,nz],[rus,rus,rus]);
if average_flag ==1
    xyz(1,:,:) = repmat(xyz(1,:,end)',1,ny);
end


% Case info
% Undersampled
Lx = xyz(1,end,1,1) - xyz(1,1,1,1);
Lz = xyz(3,1,1,end) - xyz(3,1,1,1);
x = squeeze(xyz(1,:,1,1));
y = squeeze(xyz(2,1,:,1));
z = squeeze(xyz(3,1,1,:));
dx = x(2) - x(1);
dz = z(2) - z(1);


%% Read data And undersample at the same time
% Read UVW
UVW = h5read(fname_in,'/velocity',[2,1,1], [nx,ny,nz],[rus,rus,rus]);
U = squeeze(UVW(1,:,:,:));
V = squeeze(UVW(2,:,:,:));
W = squeeze(UVW(3,:,:,:));

% Read P
P = h5read(fname_in,'/pressure',[2,1,1], [nx,ny,nz],[rus,rus,rus]);

% Read rho
rho = h5read(fname_in,'/rho',[2,1,1], [nx,ny,nz],[rus,rus,rus]);

% Read T
T = h5read(fname_in,'/temperature',[2,1,1], [nx,ny,nz],[rus,rus,rus]);



%% Filter to create filtered fields

% filter U
[U1,V1,W1] = filtervels(U,V,W,nx,ny,nz,number_filters_us);
[Uf,Vf,Wf] = filtervels(U1,V1,W1,nx,ny,nz,number_filters);
disp('filtered vels');

% filter P
% P1 = filterprodvel(P,1,nx,ny,nz,number_filters_us);
% Pf= filterprodvel(P1,1,nx,ny,nz,number_filters);
P1 = filterscalar(P,nx,ny,nz,number_filters_us,1,2,3);
Pf = filterscalar(P1,nx,ny,nz,number_filters,1,2,3);
disp('filtered P')

% filter rho
% rho1 = filterprodvel(rho,1,nx,ny,nz,number_filters_us);
% rhof= filterprodvel(rho1,1,nx,ny,nz,number_filters);
rho1 = filterscalar(rho,nx,ny,nz,number_filters_us,1,2,3);
rhof = filterscalar(rho1,nx,ny,nz,number_filters,1,2,3);
disp('filtered rho')

% filter T
% T1 = filterprodvel(T,1,nx,ny,nz,number_filters_us);
% Tf= filterprodvel(T1,1,nx,ny,nz,number_filters);
T1 = filterscalar(T,nx,ny,nz,number_filters_us,1,2,3);
Tf = filterscalar(T1,nx,ny,nz,number_filters,1,2,3);
disp('filtered T')

%% Derivatives
% duidxj
% duidxj(1,1) = dudx
% duidxj(1,2) = dudy
% duidxj(1,3) = dudz
duidxj = zeros(nx,ny,nz,3,3);
duidxj(:,:,:,1,1) = FirstDivX(Uf,nx,ny,nz,dx,0);
duidxj(:,:,:,2,1) = FirstDivX(Vf,nx,ny,nz,dx,0);
duidxj(:,:,:,3,1) = FirstDivX(Wf,nx,ny,nz,dx,0);

duidxj(:,:,:,1,2) = FirstDivY(Uf,nx,ny,nz,y);
duidxj(:,:,:,2,2) = FirstDivY(Vf,nx,ny,nz,y);
duidxj(:,:,:,3,2) = FirstDivY(Wf,nx,ny,nz,y);

duidxj(:,:,:,1,3) = FirstDivZ(Uf,nx,ny,nz,dz);
duidxj(:,:,:,2,3) = FirstDivZ(Vf,nx,ny,nz,dz);
duidxj(:,:,:,3,3) = FirstDivZ(Wf,nx,ny,nz,dz);
disp('Computed first derivatives')

% d2ui/dxjdxk
% (:,:,:,:,1) = d2uidx1dx1
% (:,:,:,:,2) = d2uidx2dx2
% (:,:,:,:,3) = d2uidx3dx3
% (:,:,:,:,4) = d2uidx2dx1
% (:,:,:,:,5) = d2uidx3dx1
% (:,:,:,:,6) = d2uidx2dx3
d2uidxjdxk = zeros(nx,ny,nz,3,6);
d2uidxjdxk(:,:,:,1,1) = FirstDivX(duidxj(:,:,:,1,1),nx,ny,nz,dx,0); %d2udxdx
d2uidxjdxk(:,:,:,2,1) = FirstDivX(duidxj(:,:,:,2,1),nx,ny,nz,dx,0); %d2vdxdx
d2uidxjdxk(:,:,:,3,1) = FirstDivX(duidxj(:,:,:,3,1),nx,ny,nz,dx,0); %d2wdxdx

d2uidxjdxk(:,:,:,1,2) = FirstDivY(duidxj(:,:,:,1,2),nx,ny,nz,y); %d2udydy
d2uidxjdxk(:,:,:,2,2) = FirstDivY(duidxj(:,:,:,2,2),nx,ny,nz,y); %d2vdydy
d2uidxjdxk(:,:,:,3,2) = FirstDivY(duidxj(:,:,:,3,2),nx,ny,nz,y); %d2wdydy

d2uidxjdxk(:,:,:,1,3) = FirstDivZ(duidxj(:,:,:,1,3),nx,ny,nz,dz); %d2udzdz
d2uidxjdxk(:,:,:,2,3) = FirstDivZ(duidxj(:,:,:,2,3),nx,ny,nz,dz); %d2vdzdz
d2uidxjdxk(:,:,:,3,3) = FirstDivZ(duidxj(:,:,:,3,3),nx,ny,nz,dz); %d2wdzdz

d2uidxjdxk(:,:,:,1,4) = FirstDivX(duidxj(:,:,:,1,2),nx,ny,nz,dx,0); %d2udydx
d2uidxjdxk(:,:,:,2,4) = FirstDivX(duidxj(:,:,:,2,2),nx,ny,nz,dx,0); %d2vdydx
d2uidxjdxk(:,:,:,3,4) = FirstDivX(duidxj(:,:,:,3,2),nx,ny,nz,dx,0); %d2wdydx

d2uidxjdxk(:,:,:,1,5) = FirstDivX(duidxj(:,:,:,1,3),nx,ny,nz,dx,0); %d2udzdx
d2uidxjdxk(:,:,:,2,5) = FirstDivX(duidxj(:,:,:,2,3),nx,ny,nz,dx,0); %d2vdzdx
d2uidxjdxk(:,:,:,3,5) = FirstDivX(duidxj(:,:,:,3,3),nx,ny,nz,dx,0); %d2wdzdx

d2uidxjdxk(:,:,:,1,6) = FirstDivZ(duidxj(:,:,:,1,2),nx,ny,nz,dz); %d2udydz
d2uidxjdxk(:,:,:,2,6) = FirstDivZ(duidxj(:,:,:,2,2),nx,ny,nz,dz); %d2vdydz
d2uidxjdxk(:,:,:,3,6) = FirstDivZ(duidxj(:,:,:,3,2),nx,ny,nz,dz); %d2wdydz
disp('Computed second derivatives')


% dTdxi
% dTdxi(1) = dTdx
% dTdxi(2) = dTdy
% dTdxi(3) = dTdz
dTdxi = zeros(nx,ny,nz,3);
dTdxi(:,:,:,1) = FirstDivX(Tf,nx,ny,nz,dx,0);
dTdxi(:,:,:,2) = FirstDivY(Tf,nx,ny,nz,y);
dTdxi(:,:,:,3) = FirstDivZ(Tf,nx,ny,nz,dz);
disp('Computed first derivatives')

% d2T/dxjdxk
% (:,:,:,1) = d2Tdx1dx1
% (:,:,:,2) = d2Tdx2dx2
% (:,:,:,3) = d2Tdx3dx3
% (:,:,:,4) = d2Tdx2dx1
% (:,:,:,5) = d2Tdx3dx1
% (:,:,:,6) = d2Tdx2dx3
d2Tdxjdxk = zeros(nx,ny,nz,6);
d2Tdxjdxk(:,:,:,1) = FirstDivX(dTdxi(:,:,:,1),nx,ny,nz,dx,0);
d2Tdxjdxk(:,:,:,2) = FirstDivY(dTdxi(:,:,:,2),nx,ny,nz,y);
d2Tdxjdxk(:,:,:,3) = FirstDivZ(dTdxi(:,:,:,3),nx,ny,nz,dz);

d2Tdxjdxk(:,:,:,4) = FirstDivX(dTdxi(:,:,:,2),nx,ny,nz,dx,0);
d2Tdxjdxk(:,:,:,5) = FirstDivX(dTdxi(:,:,:,3),nx,ny,nz,dx,0);
d2Tdxjdxk(:,:,:,6) = FirstDivZ(dTdxi(:,:,:,2),nx,ny,nz,dz);



%% Compute derived quantities

% speed of sound
gamma = 1.4;
af = sqrt(gamma*Pf./rhof);

% mu
Tref = 1; % 273 or edge temperature
Sref = 0.40417353102690834;
muref = 0.0004606654124440234;
mu = muref.*(Tf./Tref).^(3/2) .* (Tref + Sref)./(Tf + Sref);

% Pr
Pr = 0.72;

% specific heat at constant pressure
cp = 1005; % air
% R = 8.314;
% cp = gamma/(gamma-1)*R;

% Ma_edge

% Re local
Re = rhof.*(sqrt(Uf.^2 + Vf.^2 + Wf.^2)).*squeeze(xyz(2,:,:,:))./mu;

% (T - Twall)/(Te - Tw)
% = (T/Te - Tw/Te)/(1-Tw/Te)
TmTw_TemTw = (Tf - T(1,1,1))./(1 - T(1,1,1));


%%
% list of variables to pass into ML
% T/Te
% U/Ue
% V/Ve
% W/We
% rho/rhoe
% P/Pe

% c/ce
% Mae
% Ree
% Pr

% duidxj * delta_kevin/Ue
% d2ui/dxjdxk *delta_kevin^2/Ue

% mu/muref

% (T-Tw)/(Te-Tw)
% y/delta_kevin


%% Saving 

% Save full filtered field
% create file
fileID = H5F.create(fname_out,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

% Write info parameters
h5writeatt(fname_out,'/','Source file name', fname_in);
h5writeatt(fname_out,'/','Ratio of undersampling, rus', rus);
h5writeatt(fname_out,'/','number_filters', number_filters);
h5writeatt(fname_out,'/','number_filters_us', number_filters_us);
h5writeatt(fname_out,'/','fil_size', fil_size);

% Write some data
    % coordinates
h5create(fname_out,'/centerCoordinates_x',[nx,ny,nz]);
h5write(fname_out, '/centerCoordinates_x', squeeze(xyz(1,:,:,:)));
h5create(fname_out,'/centerCoordinates_y',[nx,ny,nz]);
h5write(fname_out, '/centerCoordinates_y', squeeze(xyz(2,:,:,:)));
h5create(fname_out,'/centerCoordinates_z',[nx,ny,nz]);
h5write(fname_out, '/centerCoordinates_z', squeeze(xyz(3,:,:,:)));

    % Pressure
h5create(fname_out,'/pressure',[nx,ny,nz]);
h5write(fname_out, '/pressure', Pf);

    % rho
h5create(fname_out,'/rho',[nx,ny,nz]);
h5write(fname_out, '/rho', rhof);

    % Temperature
h5create(fname_out,'/temperature',[nx,ny,nz]);
h5write(fname_out, '/temperature', Tf);

    % Velocity
h5create(fname_out,'/velocity_x',[nx,ny,nz]);
h5write(fname_out, '/velocity_x', Uf);
h5create(fname_out,'/velocity_y',[nx,ny,nz]);
h5write(fname_out, '/velocity_y', Vf);
h5create(fname_out,'/velocity_z',[nx,ny,nz]);
h5write(fname_out, '/velocity_z', Wf);

    % Sound speed
h5create(fname_out,'/soundspeed',[nx,ny,nz]);
h5write(fname_out, '/soundspeed', af);
    
    % mu
h5create(fname_out,'/mu',[nx,ny,nz]);
h5write(fname_out, '/mu', mu);



    % Velocity Derivatives
for i = 1:3
    for j = 1:3
        dataname = sprintf('/du%idx%i',i,j);
        h5create(fname_out,dataname,[nx,ny,nz]);
        h5write(fname_out, dataname, duidxj(:,:,:,i,j));
    end
end

for i = 1:3
    for j = 1:6
        dataname = sprintf('/d2u(%i)dxidxj(%i)',i,j);
        h5create(fname_out,dataname,[nx,ny,nz]);
        h5write(fname_out, dataname, d2uidxjdxk(:,:,:,i,j));
    end
end

    % Temeprature Derivatives
for i = 1:3
    dataname = sprintf('/dTdx%i',i);
    h5create(fname_out,dataname,[nx,ny,nz]);
    h5write(fname_out, dataname, dTdxi(:,:,:,i));
end

for j = 1:6
    dataname = sprintf('/d2Tdxidxj(%i)',j);
    h5create(fname_out,dataname,[nx,ny,nz]);
    h5write(fname_out, dataname, d2Tdxjdxk(:,:,:,j));
end


% fileID = H5F.create(fname_out,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
% type_id = H5T.array_create('H5T_ARRAY',3,[nx,ny,nz]);
% H5T.insert(type_id,'xyz',12,'H5T_NATIVE_DOUBLE');
% h5create(fname_out,'/centerCooridinates',[nx,ny,nz])

disp('saved filtered field')


%% Compute outputs

% Shear stress
% DNS shear stress at the undersampled locations
mu_us = muref.*(T./Tref).^(3/2) .* (Tref + Sref)./(T + Sref);
tau_wall = ComputeTauWall(U,mu_us,y);

% Heat flux
% DNS shear stress at the undersampled locations
kappa = mu/Pr*cp;
heatflux_wall = ComputeTauWall(-T,kappa,y);

disp('Computed outputs')


% return

%% old saving Save in format ready for training
if output_traning_data == 1
    % make directory for saving the files
    tmpstr=sprintf('_randperm%i',randperm_flag);
    training_folder_name = [filtered_file_name(1:end-4),tmpstr,'_training'];
    mkdir([training_out_path, training_folder_name]);

    xind_start = int32(floor(xmin*nx));
    xind_end = int32(floor(xmax*nx));
    nxtrain = xind_end - xind_start + 1;
    
    n_train_total = nz*nxtrain;
    nbatches = int32(floor(n_train_total/mini_batch_size));

    % Loop through all the mini_batches
    jtmp = 1; %xind_start;
    ktmp = 1;
    % randomly permute data
    j_randperm = xind_start:xind_end;
    k_randperm = 1:nz;
    if randperm_flag == 1
        j_randperm = j_randperm(randperm(nxtrain));
        k_randperm = k_randperm(randperm(nz));
    end

    for i = 1:nbatches
        disp(sprintf('saving batch # %i',i))
        fname_tmp = sprintf('/batch_%0.5i.hdf',i);
        fname = [training_out_path, training_folder_name, fname_tmp];
        % create file
        fileID = H5F.create(fname,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

        % Write info parameters
        noutputs = 2;
        nvars =  1 ... % y
               + 3 ... % Uf + Vf + Wf 
               + 3 ... % Pf + rhof + Tf 
               + 9 ... % duidxj
               + 18 ...% d2uidxjdxk
               + 3 ... % dTdxi
               + 6 ... % d2Tdxjdxk
               + 4 ... % af + mu + kappa + Re
               + 1;    % TmTw_TemTw
        ninputs = npoints*nvars;
        
        h5writeatt(fname,'/','Source file name', fname_in);
        h5writeatt(fname,'/','Ratio of undersampling, rus', rus);
        h5writeatt(fname,'/','number_filters', number_filters);
        h5writeatt(fname,'/','number_filters_us', number_filters_us);
        h5writeatt(fname,'/','fil_size', fil_size);
        h5writeatt(fname,'/','batch number', i);
        h5writeatt(fname,'/','num outputs per set (placed at the start of each row)', noutputs);
        h5writeatt(fname,'/','Number of wall-normal points', npoints);
        h5writeatt(fname,'/','Fraction of x where data was started', xmin);
        h5writeatt(fname,'/','Fraction of x where data was ended', xmax);
        h5writeatt(fname,'/','number of features per point used', nvars);
        h5writeatt(fname,'/','total nx of filtered dataset', nxtrain);
        h5writeatt(fname,'/','total nz of filtered dataset', nz);
        h5writeatt(fname,'/','Was the data randomly permuted (0=no)', randperm_flag);        
        h5writeatt(fname,'/','ordering for non-permuted dataset, (x,z): (1,1), (1,2), (1,3)... (2,1), (2,2)...', 1);
        h5writeatt(fname,'/','dx', dx);
        h5writeatt(fname,'/','dz', dz);
        

        % Write some data
        h5create(fname,'/data',[mini_batch_size,ninputs+noutputs]);
        h5create(fname,'/indices',[mini_batch_size,2]);
        h5create(fname,'/xz_coords',[mini_batch_size,2]);
        ind = 0;

        batch_done = 0;
        for jind = jtmp:nxtrain
            j = j_randperm(jind);

            if jind > jtmp
                ktmp = 1;
            end
            
            for kind = ktmp:nz
                k = k_randperm(kind);

                ind = ind + 1;
                if ind > mini_batch_size
                    jtmp = jind;
                    ktmp = kind;
                    batch_done = 1;
                    break
                end
                
                row = [ tau_wall(j,k), heatflux_wall(j,k),...   
                      y(1:npoints)', Uf(j,1:npoints,k), Vf(j,1:npoints,k), ...
                      Wf(j,1:npoints,k), Pf(j,1:npoints,k), ...
                      rhof(j,1:npoints,k), Tf(j,1:npoints,k),...
                      af(j,1:npoints,k), mu(j,1:npoints,k), kappa(j,1:npoints,k),...
                      Re(j,1:npoints,k), TmTw_TemTw(j,1:npoints,k)];
                for l = 1:3
                    for m = 1:3
                        row = [row, duidxj(j,1:npoints,k,l,m)];
                    end
                end
                
                for l = 1:3
                    for m = 1:6
                        row = [row, d2uidxjdxk(j,1:npoints,k,l,m)];
                    end
                end

                for l = 1:3
                    row = [row, dTdxi(j,1:npoints,k,l)];
                end
                
                for m = 1:6
                    row = [row, d2Tdxjdxk(j,1:npoints,k,m)];
                end
                            
                if size(row) ~= ninputs+noutputs
                    disp('ERROR, row size wrong!!!!!')
                    pause
                end

                h5write(fname, '/data',row,[ind,1],[1,ninputs+noutputs]);
                h5write(fname, '/indices',[j,k],[ind,1],[1,2]);
                h5write(fname, '/xz_coords',[x(j),z(k)],[ind,1],[1,2]);
            end
            
            if batch_done == 1
                break
            end
        end

        % Close the hdf5 file
        H5F.close(fileID) 
    end

end




return


% testing
figure
hold on
vals = [];
xz = [];
for i = 1:11
    tmpname = sprintf('./batch_%0.5i.hdf',i);
    xz = [xz; h5read(tmpname,'/xz_coords')];
    vals = [vals; h5read(tmpname,'/data')];
end
[xq,yq] = meshgrid(x, z); % interpolation points grid
% vq = griddata(xz(:,1),xz(:,2),vals(:,9),xq,yq,'natural'); %U
    vq = griddata(xz(:,1),xz(:,2),vals(:,2),xq,yq,'natural');
Z = squeeze(vq);
contourf(xq,yq,Z);



