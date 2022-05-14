%----------------------------------------------%
%       compute statistics from fields         %
%----------------------------------------------%
clear all

%-----------------------------------------------------------------%
iplot = 0; % plot data
isave = 1; % save data in .mat
isymmetrize   = 1; % 0-> no symmetrize, 1->symmetrize respect to centerline
fname   = ['395_DNS/DNS_395_201'];
fileout = ['match_resolved_stress.mat']; 
vv_  = 300000:1000:329000;

%-------------------------------------------------------%
% read fields
ifirst = 1;
iNaN   = 0;
nfail  = 0;
nfils = 0; 
for i = vv_(1:1:end)

  nfils = nfils + 1 ; 
  fin = [fname,'.',num2str(i)]
%  fin = fname; 
  fid = fopen(fin,'r','b');

   if fid<0, return, end 
  disp(fin)

  n   = fread(fid,1,'int32');
  x   = fread(fid,n,'float64');

  n   = fread(fid,1,'int32');
  y   = fread(fid,n,'float64');

  n   = fread(fid,1,'int32');
  z   = fread(fid,n,'float64');

  n   = fread(fid,1,'int32');
  xm  = fread(fid,n,'float64');

  n   = fread(fid,1,'int32'); 
  ym  = fread(fid,n,'float64');

  n   = fread(fid,1,'int32');
  zm  = fread(fid,n,'float64');
  
  n   = fread(fid,3,'int32');
  U   = fread(fid,n(1)*n(2)*n(3),'float64');
  U   = reshape(U,n(1),n(2),n(3));

  n   = fread(fid,3,'int32');
  V   = fread(fid,n(1)*n(2)*n(3),'float64');
  V   = reshape(V,n(1),n(2),n(3));

  n   = fread(fid,3,'int32');
  W   = fread(fid,n(1)*n(2)*n(3),'float64');
  W   = reshape(W,n(1),n(2),n(3));		       

  n   = fread(fid,3,'int32');
  P   = fread(fid,n(1)*n(2)*n(3),'float64');
  P   = reshape(P,n(1),n(2),n(3));

  n_2    = fread(fid,3,'int32');
  nu_t = fread(fid,n(1)*n(2)*n(3),'float64');
  nu_t = squeeze(reshape(nu_t,n(1),n(2),n(3)));		       
  
  X = fread(fid,'float64');

xg(2:length(xm)+1,1) = xm;
  xg(1,1) = xm(1) - 2*(xm(1)-x(1));
  xg(length(xm)+2,1) = xm(length(xm)) + 2*(x(length(x))-xm(length(xm)));
  
  yg(2:length(ym)+1,1) = ym;
  yg(1,1) = ym(1) - 2*(ym(1)-y(1));
  yg(length(ym)+2,1) = ym(length(ym)) + 2*(y(length(y))-ym(length(ym)));
  zg(2:length(zm)+1,1) = zm;
  zg(1,1) = zm(1) - 2*(zm(1)-z(1));
  zg(length(zm)+2,1) = zm(length(zm)) + 2*(z(length(z))-zm(length(zm)));
  
  fclose(fid);
%   return;
  
  % stats
if i ==vv_(1)
[U_int_xyz,V_int_xyz,W_int_xyz,~,~] = ...
    interpolate_flow_fields(U,V,W,P,nu_t,y,ym);
else
[U2,V2,W2,~,~] = ...
    interpolate_flow_fields(U,V,W,P,nu_t,y,ym);
U_int_xyz = U_int_xyz + U2;
V_int_xyz = V_int_xyz + V2;
W_int_xyz = W_int_xyz + W2;

end

end
clear U  V W
U = U_int_xyz(2:end-1,:,2:end-1);
V = V_int_xyz(2:end-1,:,2:end-1);
W = W_int_xyz(2:end-1,:,2:end-1);

%-----------------------------------------------------------------%
% average over fields
nfils
U     = U/nfils;
V     = V/nfils;
W     = W/nfils;
Nx = 256;
Ny = 201;
Nz = 256;
x = x(2:end-1);
z = z(2:end-1);
save('Re395channelavg.mat','U','V','W','Nx','Ny','Nz','x','y','z')
