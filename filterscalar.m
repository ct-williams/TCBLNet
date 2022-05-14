function [U1] = filterscalar(U1,nx,ny,nz,nfilters,xflag,yflag,zflag)
for i = 1: nfilters
    if(i==1)
        U1 = filter_vasilyevxz(U1,nx,ny,nz, xflag);
        U1 = filter_vasilyevxz(U1,nx,ny,nz,zflag);
        U1 = filter_vasilyevy(U1,nx,ny,nz,yflag);
    else 
        U1 = filter_vasilyevxz(U1,nx,ny,nz, xflag);
        U1 = filter_vasilyevxz(U1,nx,ny,nz,zflag);
        U1 = filter_vasilyevy(U1,nx,ny,nz,yflag);
    end
end
end