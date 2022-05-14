function [U,V,W] = filtervels(U,V,W,nx,ny,nz,nfilters)
   for i = 1 : nfilters
       if(i==1)
            U = filter_vasilyevxz(U,nx,ny,nz, 1);
            U = filter_vasilyevxz(U,nx,ny,nz,3);
            U = filter_vasilyevy(U,nx,ny,nz,1);

            V = filter_vasilyevxz(V,nx,ny,nz, 1);
            V = filter_vasilyevxz(V,nx,ny,nz,3);
            V = filter_vasilyevy(V,nx,ny,nz,2);

            W = filter_vasilyevxz(W,nx,ny,nz, 1);
            W = filter_vasilyevxz(W,nx,ny,nz,3);
            W = filter_vasilyevy(W,nx,ny,nz,3);
          %  U=U3;
          %  V=V3;
          %  W=W3;


        else 
            U = filter_vasilyevxz(U,nx,ny,nz, 1);
            U = filter_vasilyevxz(U,nx,ny,nz,3);
            U = filter_vasilyevy(U,nx,ny,nz,1);

            V = filter_vasilyevxz(V,nx,ny,nz, 1);
            V = filter_vasilyevxz(V,nx,ny,nz,3);
            V = filter_vasilyevy(V,nx,ny,nz,2);

            W = filter_vasilyevxz(W,nx,ny,nz, 1);
            W = filter_vasilyevxz(W,nx,ny,nz,3);
            W = filter_vasilyevy(W,nx,ny,nz,3);
       end
   end
end

