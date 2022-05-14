function [U_int_xyz,V_int_xyz,W_int_xyz] = interpolate_flow_fields(U,V,W,y,ym)
    yg(2:length(ym)+1,1) = ym;
    yg(1,1) = ym(1) - 2*(ym(1)-y(1));
    yg(length(ym)+2,1) = ym(length(ym)) + 2*(y(length(y))-ym(length(ym)));

    weight_0  = (yg(2:end) - y)./(yg(2:end) - yg(1:end-1));
    weight_1  = 1 - weight_0;

    % Interpolating in y  
    for j = 1:length(y)
        U_int_y(:,j,:) = weight_0(j)*U(:,j,:)+weight_1(j)*U(:,j+1,:);
        W_int_y(:,j,:) = weight_0(j)*W(:,j,:)+weight_1(j)*W(:,j+1,:);
 %       T_int_y(:,j,:) = weight_0(j)*T(:,j,:)+weight_1(j)*T(:,j+1,:);
 %       P_int_y(:,j,:) = weight_0(j)*P(:,j,:)+weight_1(j)*P(:,j+1,:);
 %       kappa_t_int_y(:,j,:) = weight_0(j)*kappa_t_int_y(:,j,:)+weight_1(j)*kappa_t_int_y(:,j+1,:);
 %       nu_t_int_y(:,j,:) = weight_0(j)*nu_t(:,j,:)+weight_1(j)*nu_t(:,j+1,:);
    end
    V_int_y = V;
    
    % Interpolating in z
    U_int_yz = (U_int_y(:,:,2:end)+U_int_y(:,:,1:end-1))/2;
    V_int_yz = (V_int_y(:,:,2:end)+V_int_y(:,:,1:end-1))/2;
 %   T_int_yz = (T_int_y(:,:,2:end)+T_int_y(:,:,1:end-1))/2;
 %   P_int_yz = (P_int_y(:,:,2:end)+P_int_y(:,:,1:end-1))/2;
 %   kappa_t_int_yz = (kappa_t_int_y(:,:,2:end)+kappa_t_int_y(:,:,1:end-1))/2;
 %   nu_t_int_yz = (nu_t_int_y(:,:,2:end)+nu_t_int_y(:,:,1:end-1))/2;
    W_int_yz = W_int_y;
    
    % Interpolating in x
    U_int_xyz = U_int_yz;
    V_int_xyz = (V_int_yz(2:end,:,:)+V_int_yz(1:end-1,:,:))/2;
 %   T_int_xyz = (T_int_yz(2:end,:,:)+T_int_yz(1:end-1,:,:))/2;
 %   P_int_xyz = (P_int_yz(2:end,:,:)+P_int_yz(1:end-1,:,:))/2;
 %   kappa_t_int_xyz = (kappa_t_int_yz(2:end,:,:)+kappa_t_int_yz(1:end-1,:,:))/2;
 %   nu_t_int_xyz = (nu_t_int_yz(2:end,:,:)+nu_t_int_yz(1:end-1,:,:))/2;
    W_int_xyz = (W_int_yz(2:end,:,:)+W_int_yz(1:end-1,:,:))/2;
end