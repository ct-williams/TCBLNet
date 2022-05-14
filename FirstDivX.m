function der = FirstDivX(field,Nx,ny,Nz,dx,periodic_flag)
    % First order derivative | 2nd order central difference
    % Field - velocity field U, V, or W
%    field = zeros(Nx+2,ny,Nz);
%    field(1,:,:) = field(end,:,:);
%    field(end,:,:) = field(1,:,:);
%    field(2:end-1,:,:) = field;

    der = zeros(Nx,ny,Nz);

    % interior
    der = [(field(3:end,:,:)-field(1:end-2,:,:))/ (2*dx)] ;
        
    if periodic_flag == 1
        der = [    (field(2,:,:)-field(end,:,:))/ (2*dx) ; ...
                der; ...
                (field(1,:,:)-field(end-1,:,:))/(2*dx) ] ;
   else
        der = [    (field(2,:,:)-field(1,:,:))/ (dx) ; ...
                der; ...
                (field(end,:,:)-field(end-1,:,:))/(dx) ] ;
    end
 end
