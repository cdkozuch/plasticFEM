clear
close all

nten = 6; %number of unique stress/strain tensor values

%time values
dt = 0.05; tmin = 0; tmax = 5;
tvec = tmin:dt:tmax;
nt = numel(tvec);

%elastic material properties
E = 15; nu = 0.25;
lambda = (E*nu)/((1+nu)*(1-(2*nu))); %first lame constant
mu = E/(2*(1+nu)); %second lame constant

%plastic material properties
Hp = 4.5;
K0 = 35;
Kp = 2.5;

matprops = [lambda,mu,K0,Kp,Hp];

%strain evolution function
eps1 = 10; eps2 = 5;

%preallocate
ep = zeros(nten,nt); %plastic strain
alpha = zeros(1,nt); %equivalent plastic strain
qbar = zeros(nten,nt); %back stress
sigma = zeros(nten,nt); %stress
C = zeros(nten,nten,nt); %material stiffness
epsilon = zeros(nten,nt); %strain
C(:,:,1) = matStiffTen3D(lambda,mu);

%perform return mapping for each time step
for i=2:nt
    t = tvec(i);
    epsilon(:,i) = eps1*t*[1 1 3/5 0 0 0]' + eps2*sin(t)*[0 1 5/3 4 0 0]';
    [ep(:,i),sigma(:,i),qbar(:,i),alpha(i),C(:,:,i)] = ...
        returnMap(epsilon(:,i),ep(:,i-1),qbar(:,i-1),alpha(i-1),matprops);
end

%plot stress
figure
plot(tvec,sigma(1,:),tvec,sigma(2,:),...
     tvec,sigma(3,:),tvec,sigma(4,:),'LineWidth',2)
legend('s11','s22','s33','s12')
xlabel('time')
ylabel('sigma')

%plot alpha
figure
plot(tvec,alpha,'LineWidth',2)
xlabel('time')
ylabel('alpha')

%plot qbar
figure
plot(tvec,qbar(1,:),tvec,qbar(2,:),...
     tvec,qbar(3,:),tvec,qbar(4,:),'LineWidth',2)
legend('q11','q22','q33','q12') %using equation 4.3.1a
xlabel('time')
ylabel('qbar')

%plot C
figure
plot(tvec,reshape(C(1,1,:),1,nt),tvec,reshape(C(2,2,:),1,nt),...
     tvec,reshape(C(4,4,:),1,nt),'LineWidth',2)
legend('C1111','C2222','C1212') %using equation 4.3.1b
xlabel('time')
ylabel('C')
