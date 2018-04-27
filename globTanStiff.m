function [KT,Fint,epsilon2,ep2,sigma2,qbar2,alpha2,D2] = globTanStiff(elements,coors,U,ep1,qbar1,alpha1,matprops)
%%GLOBTANSTIFF computes global tangent stiffness and residual force
%Inputs:
%   elements (nel x nnel): element connectivity
%   coors (nnode x ndofn): nodal coordinates
%   U (nten x 1): displacement
%   ep1 (nten x 1): current plastic strain
%   qbar1 (nten x 1): current back stress
%   alpha1 (scalar): current equivalent plastic strain
%   matprops (1 x 5): material properties [lambda,mu,K0,Kp,Hp]
%     lambda: first lame parameter
%     mu: second lame parameter
%     K0: yield stress
%     Kp: isotropic hardening modulus
%     Hp: kinematic hardening modulus
%Ouputs:
%   KT (ndoft x ndoft): global tangent stiffness matrix
%   Fint (ndoft x 1): internal force vector
%   epsilon2 (nten x nintpt): strain
%   ep2 (nten x nintpt): new plastic strain
%   sigma2 (nten x nintpt): new stress
%   qbar2 (nten x nintpt): new back stress
%   alpha2 (1 x nintpt): new equivalent plastic strain
%   D2 (nten x nten x nintpt): new material stiffness tensor

nten = 6; %number of unique stress/strain values in 3D

ndofn = size(coors,2); %dofs per node
nnode = size(coors,1); %number of nodes
nnel = size(elements,2); %nodes per elements
nel = size(elements,1); %number of elements
ndoft = nnode*ndofn; %total number of degrees of freedom
ndofel = nnel*ndofn; %number of degrees of freedom per element

gp = (1/sqrt(3))*[-1 1]; %gauss points
ngp = numel(gp); %number of guass points
nintpel = ngp^ndofn; %number of integration points per element
nintpt = nintpel*nel; %number of total integration points

%preallocate return variables
epsilon2 = zeros(nten,nintpt);
ep2 = zeros(nten,nintpt);
sigma2 = zeros(nten,nintpt);
qbar2 = zeros(nten,nintpt);
alpha2 = zeros(1,nintpt);
D2 = zeros(nten,nten,nintpt);

KT = spalloc(ndoft,ndoft,ndoft); %global stiffness matrix
Fint = spalloc(ndoft,1,ndoft); %residual force vector

%create global stiffness and tangent stiffness matrices
intpi = 1; %indegration point index
for iter=1:nel
    ind = globElemInd(elements,iter,ndofn); %global indices for element
    ke = zeros(ndofel); %element stiffness matrix
    fe = zeros(ndofel,1); %element internal force vector
    
    for i=1:ngp
      xi = gp(i); %first natural coordinate
      
      for j=1:ngp
        eta = gp(j); %second natural coordinate
        
        for k=1:ngp
          zeta = gp(k); %third natural coordinate

          [B,J] = bmatHex8(coors,[xi,eta,zeta]); %B matrix   
          epsilon2(:,intpi) = B*U(ind);
          
          [ep2(:,intpi),sigma2(:,intpi),qbar2(:,intpi),alpha2(intpi),...
            D2(:,:,intpi)] = returnMap(epsilon2(:,intpi),ep1(:,intpi),...
            qbar1(:,intpi),alpha1(intpi),matprops); %return mapping
           
          detJ = det(J);
          ke = ke + detJ*B'*D2(:,:,intpi)*B;
          fe = fe + detJ*B'*sigma2(:,intpi);
          intpi = intpi + 1;
        end
      end
    end

    Fint(ind) = Fint(ind) + fe; %update global residual force vector
    KT(ind,ind) = KT(ind,ind) + ke; %update global stiffness matrix
end

end