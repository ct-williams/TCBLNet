function [U1] = filtertensor(U1,nx,ny,nz,nfilters)
for i = 1: nfilters
    if(i==1)
        U1 = filter_vasilyevxz(U1,nx,ny,nz, 1);
        U1 = filter_vasilyevxz(U1,nx,ny,nz,3);
        U1 = filter_vasilyevy(U1,nx,ny,nz,1);
        else 
        U1 = filter_vasilyevxz(U1,nx,ny,nz, 1);
        U1 = filter_vasilyevxz(U1,nx,ny,nz,3);
        U1 = filter_vasilyevy(U1,nx,ny,nz,1);
         end
end
end
