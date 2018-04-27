function [B,J] = bmatHex8(coors,natcoors)
%%BMATHEX8 computes strain-displacement matrix for 8-node hex
%Input:
%   coors (8x3 matrix): coordinates of nodes
%   natcoors (1x3 vector): natural coordinates [xi,eta,zeta]
%Ouput:
%   B (6x24 matrix): strain-displacement matrix
%   J (3x3 matrix): Jacobian

nnel = size(coors,1); %number of nodes per element
ndofn = size(coors,2); %number of degrees of freedom per node
xi = natcoors(1); eta = natcoors(2); zeta = natcoors(3);

%Calculation of shape funciton derivatives
dNdr = 0.125*[ -1*(1-eta)*(1-zeta), (1-eta)*(1-zeta), ...
               (1+eta)*(1-zeta), -1*(1+eta)*(1-zeta), ...
               -1*(1-eta)*(1+zeta), (1-eta)*(1+zeta), ...
               (1+eta)*(1+zeta), -1*(1+eta)*(1+zeta) ];                    
dNds = 0.125*[ -1*(1-xi)*(1-zeta), -1*(1+xi)*(1-zeta), ...
               (1+xi)*(1-zeta), (1-xi)*(1-zeta), ...
               -1*(1-xi)*(1+zeta), -1*(1+xi)*(1+zeta), ...
               (1+xi)*(1+zeta), (1-xi)*(1+zeta) ];
dNdt = 0.125*[ -1*(1-xi)*(1-eta), -1*(1+xi)*(1-eta), ...
               -1*(1+xi)*(1+eta), -1*(1-xi)*(1+eta), ...
               (1-xi)*(1-eta), (1+xi)*(1-eta), ...
               (1+xi)*(1+eta), (1-xi)*(1+eta) ];
dN = [dNdr;dNds;dNdt]; %shape function derivatives

%Calculation of Jacobian
J = zeros(ndofn);
for iter=1:nnel
    icoors = [coors(iter,1) coors(iter,2) coors(iter,3)];
    J = J + [ dN(1,iter)*icoors;
              dN(2,iter)*icoors;
              dN(3,iter)*icoors ];
end
Jinv = inv(J); %inverse of Jacobian

%Calculation of displacement-strain matrix
for iter=1:nnel
    DNi = Jinv*dN(:,iter); 
    ind = (iter-1)*3 + 1;
    B(:,ind:ind+2) = [ DNi(1) 0 0; 
                       0 DNi(2) 0; 
                       0 0 DNi(3);
                       DNi(2) DNi(1) 0;
                       0 DNi(3) DNi(2);
                       DNi(3) 0 DNi(1) ];
end

end