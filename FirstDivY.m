function dfdy = FirstDivY(field,Nx,ny,Nz,y)
    % First order derivative | 2nd Order accurate due to the mesh being
    % analytically stretched
    % Field - velocity field U, V, or W
    dfdy = zeros(Nx,ny,Nz);
    for i = 2 : ny-1 
    dfdy(:,i,:) = (field(:,i+1,:)-field(:,i-1,:))./(y(i+1)-y(i-1)); 
    end   
    dfdy(:,1,:) = (field(:,2,:)-field(:,1,:))/(y(2)-y(1)); % One sided differences (1st order)
    dfdy(:,end,:) = (field(:,end,:)-field(:,end-1,:))/(y(end)-y(end-1)); % One sided differences (1st order)
end
