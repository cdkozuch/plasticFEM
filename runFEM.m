%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script executes the non-linear finite element method for a mesh 
% composed of one 8-node hexahedral element with combined isotropic and
% kinematic hardening. The element is placed under uni-axial tension by 
% applying a fixed dispalacement corresponding to 2% strain in the e1 
% direction to nodes 2,3,6,7 while pinning node 1, allowing only e2
% displacement in node 4, and preventing only e1 displacement in nodes 5,8.
%
%      5-------8
%    / |     / |
%   6--|----7  |
%   |  1----|--4
%   | /     | /
%   2-------3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% Get nodal coordinates from data files
load 'nodes.dat';
load 'elements.dat';
nnode = size(nodes,1); %number of nodes 
nel = size(elements,1); %number of elements
nnel = size(elements,2); %number of nodes per element
ndofn = size(nodes,2); %number of degrees of freedom per node
coors = zeros(nnel,ndofn,nel); %nodal coordinates for each element
ndoft = nnode*ndofn; %total number of degrees of freedom
ndofel = nnel*ndofn; %number of degrees of freedom per element

for i=1:nel
  for j=1:nnel
    coors(j,:,i) = nodes(elements(i,j),:);
  end
end

%% elastic material properties
E = 11000; %elastic modulus
nu = 0.25; %poisson's ratio
lambda = (E*nu)/((1+nu)*(1-(2*nu))); %first lame constant
mu = E/(2*(1+nu)); %second lame constant

%% Plastic material properties
K0 = 120;
Kp = 900;
Hp = 2.5;
matprops = [lambda,mu,K0,Kp,Hp];
      
%% Define boundary conditions
utotal = 0.02; %total displacement of top surface in loading direction
nstep1 = 100; %number of steps until displacement total
nstep2 = 120; %final number of steps
du = utotal/nstep1; %displacement increment
fixDispDofs = [1 2 3 10 12 13 22 4 7 16 19]; %fixed displacement DoFs
nFixDispDofs = numel(fixDispDofs);
dFixDispVals = [zeros(1,7) du*ones(1,4)]; %fixed displacment increment values
freeDispDofs = setdiff(1:ndoft,fixDispDofs);
Fext = spalloc(ndoft,1,ndoft); %no external forces

%% Initialize plastic history variables
nintpel = 8; %number of integration points per element (2x2x2 quadrature)
nintpt = nintpel*nel; %number of integration points
nten = 6; %number of unique strain/stress tensor elements
ep = zeros(nten,nintpt); %plastic strain
qbar = zeros(nten,nintpt); %kinematic hardening
alpha = zeros(nintpt); %effective plastic strain

%% Preallocate history storage
sigmaRec = cell(1); %stress history
sigmaRec{1} = zeros(nten,nintpt);
epsilonRec = cell(1); %strain history
epsilonRec{1} = zeros(nten,nintpt);
epRec = cell(1); %plastic strain history
epRec{1} = ep;
alphaRec = zeros(nstep2+1,nintpt); %effective plastic strain history
qbarRec = cell(1); %kinematic hardening history
qbarRec{1} = qbar;
Urec = zeros(ndoft,nstep2+1); %displacement (columns are time steps)    
Crec = cell(1); %material stiffness tensor history
C = matStiffTen3D(lambda,mu);
Crec{1} = repmat(C,1,1,nintpt);

%% Step through analysis
maxni = 200; %maximum number of newton iterations
tolFac = 1e-12; %tolerance factor

tic
for t=1:nstep2
    U = zeros(ndoft,1);
    U(fixDispDofs) = t*dFixDispVals;
    
    %Newton iterations
    for ni=1:maxni
      %Build tangent stiffness and internal force
      [KT,Fint,epsilon,ep,sigma,qbar,alpha,C] = ...
        globTanStiff(elements,coors,U,epRec{t},qbarRec{t},alphaRec(t,:),matprops);     
      
      %Compute residual force
      RF = Fext - Fint;
      
      %Apply boundary conditions
      KT(fixDispDofs,:) = zeros(nFixDispDofs,ndoft);
      KT(fixDispDofs,fixDispDofs) = eye(nFixDispDofs);

      %Evaluate condition
      res = norm(RF(freeDispDofs));
      fprintf('Time step: %i;\tIteration: %i;\tResidual: %g\n',t,ni,res)   
      if ni==1
          tol = tolFac*res;
      elseif res <= tol
        break
      elseif ni==maxni
        error('Convergence failure!')
      end

      %Update displacement guess
      U(freeDispDofs) = U(freeDispDofs) + ...
        KT(freeDispDofs,freeDispDofs)\RF(freeDispDofs);
    end 

    %record history variables
    Urec(:,t+1) = U;
    epRec{t+1} = ep;
    epsilonRec{t+1} = epsilon;
    sigmaRec{t+1} = sigma;
    Crec{t+1} = C;
    qbarRec{t+1} = qbar;
    alphaRec(t+1,:) = alpha;
end

%% Plot results
tf1 = 1; dt = tf1/nstep1; t0 = 0; tf2 = dt*nstep2;
tvec = t0:dt:tf2; %time vector
intpi = 2; %integration point at which to plot

%plot stress
s11 = zeros(1,nstep2); s22 = zeros(1,nstep2);
s33 = zeros(1,nstep2); s12 = zeros(1,nstep2);
for t=1:nstep2+1
  s11(t) = sigmaRec{t}(1,intpi); s22(t) = sigmaRec{t}(2,intpi);
  s33(t) = sigmaRec{t}(3,intpi); s12(t) = sigmaRec{t}(4,intpi);
end
figure
plot(tvec,s11,tvec,s22,tvec,s33,tvec,s12,'LineWidth',2)
legend('s11','s22','s33','s12','Location','NorthWest')
xlabel('time')
ylabel('sigma')

%plot alpha
figure
plot(tvec,alphaRec(:,intpi),'LineWidth',2)
xlabel('time')
ylabel('alpha')

%plot qbar
q11 = zeros(1,nstep2); q22 = zeros(1,nstep2); 
q33 = zeros(1,nstep2); q12 = zeros(1,nstep2);
for t=1:nstep2+1
  q11(t) = qbarRec{t}(1,intpi); q22(t) = qbarRec{t}(2,intpi);
  q33(t) = qbarRec{t}(3,intpi); q12(t) = qbarRec{t}(4,intpi);
end
figure
plot(tvec,q11,tvec,q22,tvec,q33,tvec,q12,'LineWidth',2)
legend('qbar11','qbar22','qbar33','qbar12','Location','NorthWest')
xlabel('time')
ylabel('qbar')

%plot C
C1111 = zeros(1,nstep2); C2222 = zeros(1,nstep2); C1212 = zeros(1,nstep2);
for t=1:nstep2+1
  C1111(t) = Crec{t}(1,1,intpi);
  C2222(t) = Crec{t}(2,2,intpi);
  C1212(t) = Crec{t}(4,4,intpi);
end
figure
plot(tvec,C1111,tvec,C2222,tvec,C1212,'LineWidth',2)
legend('C1111','C2222','C1212')
xlabel('time')
ylabel('C')


