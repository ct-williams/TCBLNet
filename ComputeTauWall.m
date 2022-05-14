function [tau_wall] = ComputeTauWall(U,mu,y)
    % Computes the wall shear stress at the wall
    % (just at j=1, Note that HTR has a value on the 
    % wall)
    % (I am using a first order scheme for now)

    dudy_wall = (U(:,2,:) - U(:,1,:))./(y(2)-y(1)); % First order
    tau_wall = mu(:,1,:).*dudy_wall;
    tau_wall = squeeze(tau_wall);
end
