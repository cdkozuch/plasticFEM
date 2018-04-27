function [ep2,sigma2,qbar2,alpha2,C2] = returnMap(epsilon2,ep1,qbar1,alpha1,matprops)
%%RETURNMAP performs return mapping routine
%note: nten = 6 for 3D
%Input:
%   epsilon2 (nten x 1): new strain
%   ep1 (nten x 1): current plastic strain
%   qbar1 (nten x 1): current back stress
%   alpha1 (scalar): current equivalent plastic strain
%   matprops (1 x 5): material properties [lambda,mu,K0,Kp,Hp]
%     lambda: first lame parameter
%     mu: second lame parameter
%     K0: yield stress
%     Kp: isotropic hardening modulus
%     Hp: kinematic hardening modulus
%Output:
%   ep2 (nten x 1): new plastic strain
%   sigma2 (nten x 1): new stress
%   qbar2 (nten x 1): new back stress
%   alpha2 (scalar): new equivalent plastic strain
%   C2 (nten x nten): new material stiffness matrix 

lambda = matprops(1); mu = matprops(2); kappa = lambda + (2/3)*mu;
K0 = matprops(3); Kp = matprops(4); Hp = matprops(5);

treps = sum(epsilon2(1:3)); %trace of epsilon
vecIden2 = [ones(3,1);zeros(3,1);]; %vector form of 2-tensor identity matrix (eq. 4.3.1c)
hydeps = (1/3)*treps*vecIden2; %hydrostatic part of epsilon
deveps = epsilon2 - hydeps; %deviatoric part of epsilon

strial = 2*mu*(deveps - ep1); %trial deviatoric stress
xitrial = strial - qbar1; %relative stress 
xitrimag = sqrt(sum(xitrial(1:3).^2) + sum(2*(xitrial(4:6).^2)));

ftrial = xitrimag - sqrt(2/3)*(K0 + Kp*alpha1); %yield function
if ftrial <= 0
    %still in elastic zone
    ep2 = ep1; sigma2 = 3*kappa*hydeps + strial;
    qbar2 = qbar1; alpha2 = alpha1;
    C2 = matStiffTen3D(lambda,mu);  
    return
end

%update alpha
dgamma = ftrial/(2*mu + (2/3)*(Kp + Hp)); %consistency parameter
alpha2 = alpha1 + sqrt(2/3)*dgamma;
 
%update back-stress
nvec = xitrial/xitrimag; %directional vector
qbar2 = qbar1 + (2/3)*Hp*dgamma*nvec;

%update plastic strain
ep2 = ep1 + dgamma*nvec;

%compute stress
sigma2 = 3*kappa*hydeps + strial - 2*mu*dgamma*nvec;

%compute stiffness tensor
theta = 1 - 2*mu*dgamma/xitrimag;
thetabar = 1/(1 + (Kp + Hp)/(3*mu)) - (1 - theta);
matIden4 = diag([ones(1,3),0.5*ones(1,3)]); %4-tensor identity in 2-tensor form (eq. 4.3.1c)
II = [ones(3),zeros(3);zeros(3),zeros(3)]; %dyadic product of two 2-tensor identity
C2 = kappa*II + 2*mu*theta*(matIden4 - (1/3)*II) - 2*mu*thetabar*(nvec*nvec');

end


