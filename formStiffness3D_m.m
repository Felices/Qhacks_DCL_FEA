function [K_element] = formStiffness3D_m(local_nodecoord, D)

% Blank Matrix
K_element = zeros(24);
% Gauss Points
GaussPoint = [-0.577350269189626, 0.577350269189626];
% Because I will only ever use 2 point Gaussian Quadrature and the weights
% are all one, I didn't include the weights

elenodecoordx = local_nodecoord(:, 1);
elenodecoordy = local_nodecoord(:, 2);
elenodecoordz = local_nodecoord(:, 3);
for p = 1:length(GaussPoint)
    for q = 1:length(GaussPoint)
        for r = 1:length(GaussPoint)
            xi = GaussPoint(p);
            eta = GaussPoint(q);
            sigma = GaussPoint(r);
            % Shape Function derivatives
            dShape = (1/8) * [  -(1-eta)*(1-sigma), -(1-eta)*(1+sigma), -(1+eta)*(1+sigma), -(1+eta)*(1-sigma),   (1-eta)*(1-sigma),  (1-eta)*(1+sigma), (1+eta)*(1+sigma),  (1+eta)*(1-sigma);
                                -(1-xi)*(1-sigma),  -(1-xi)*(1+sigma),   (1-xi)*(1+sigma),   (1-xi)*(1-sigma),   -(1+xi)*(1-sigma),  -(1+xi)*(1+sigma),  (1+xi)*(1+sigma),   (1+xi)*(1-sigma);
                                -(1-xi)*(1-eta),     (1-xi)*(1-eta),     (1-xi)*(1+eta),    -(1-xi)*(1+eta),     -(1+xi)*(1-eta),     (1+xi)*(1-eta),    (1+xi)*(1+eta),    -(1+xi)*(1+eta) ]; 
            % Jacobians
            Jacob = [dShape(1,:) * elenodecoordx, dShape(1,:) * elenodecoordy, dShape(1,:) * elenodecoordz;
                     dShape(2,:) * elenodecoordx, dShape(2,:) * elenodecoordy, dShape(2,:) * elenodecoordz;
                     dShape(3,:) * elenodecoordx, dShape(3,:) * elenodecoordy, dShape(3,:) * elenodecoordz];
            invJacobian = inv(Jacob);
            disp(invJacobian);
            disp(dShape);
            XYZderivatives =  invJacobian * dShape;
            
            % B Matrix
            B = [XYZderivatives(1,:), zeros(1,8),          zeros(1,8);
                 zeros(1,8),          XYZderivatives(2,:), zeros(1,8);
                 zeros(1,8),          zeros(1,8),          XYZderivatives(3,:);
                 XYZderivatives(2,:), XYZderivatives(1,:), zeros(1,8);
                 zeros(1,8),          XYZderivatives(3,:), XYZderivatives(2,:);
                 XYZderivatives(3,:), zeros(1,8),          XYZderivatives(1,:)];
             
            % Calculate the Stiffness Matrix
            K_element = K_element +( B' * D * B * det(Jacob));
        end
    end
end
end




