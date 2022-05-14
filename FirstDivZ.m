function der = FirstDivZ(field,Nx,ny,Nz,dz)
    % First order derivative | 2nd Order central difference
    % Field - velocity field U, V, or W
    %field = zeros(Nx,ny,Nz+2);
    %field(:,:,1) = field(:,:,end);
    %field(:,:,end) = field(:,:,1);
    %field(:,:,2:end-1) = field;
    %dfdz = zeros(Nx,ny,Nz);
%    dfdz = (field(:,:,3:end)-field(:,:,1:end-2))/dz;


      der = [ permute(  (field(:,:,2)-field(:,:,end))/ (2*dz) , [3 1 2] ) ; ...
                permute( (field(:,:,3:end)-field(:,:,1:end-2))/ (2*dz) , [ 3 1 2 ] );  ...
                permute( (field(:,:,1)-field(:,:,end-1))/(2*dz) , [ 3 1 2] ) ] ;
      der = permute(der, [2 3  1 ] ) ; 


end
