clc;clear all;close all;
%{
 A*Nn = B*ft where 
 -  A is a 6x6 matrix (could be 6x12 if the second order terms are added)
    A is the calibration matrix from the transducer signal to the forces 
    at a cross-section (shaftObj.csF)

 -  Nn is a 6x1 matrix (could be 12x1 if the second order terms are added)
    Nn is the vecotr of normalized difference over sum signals.

 -  B is the mapping from the tip forces (shaftObj.tpF) to the forces at
    a cross-section of the shaft with distance Lp from its clamping point.
    The structure of B matrix is derived in bendingModel.pptx 

    B is dependent on the the instrument insertion into the trocar (Ls),
    the location of the cross section (Lp), the length of the shaft
    (shaftObj.L), and the physical properties of the shaft lumped into a
    single parameter c = 6*E*I/ks.

    During calibration of the sensor when installed on an instrument, the
    known parameters are Nn, ft, and Ls for every measurement. The goal of
    calibration is to find the A matrix, L, Lp and ks. 

    when A is a 6x6 matrix, the vector of optimization parameters (x) is
    40x1 in which the first 36 elements are the elements of the A matrix
    x(1 ... 36) --> A (6x6)
    x(37)       --> L
    x(38)       --> cx
    x(39)       --> cy
    x(40)       --> Lp
    x(41)       --> offsetLs
%}
%% Initialize the shaft properties
E = 100E9;            % modulus of elasticity (N/m2)
G = 3E9;              % shear modulus of elasticity (N/m2)

ri = 5.99/2/1000;     % inner radius of the shaft (m)
ro = 8.36/2/1000;     % outer radius of the shaft (m)
L = (430)/1000;       % length of the shaft (m)
Lp = 50e-3;           % Lp is the distance from the base of the instrument
                      % at which the sensor is mounted in (m)
%% ks is the spring stiffness in the trocar
ks = 10000;             % ks is the trocar stiffness in N/m
%% construct a shaft object
shaftObj = shaft(L,ri,ro,E,G,ks);
%% read the force data applied at the tip of the instrument
ft = readForceData();
ft(:,4:end) = ft(:,4:end)/1000;  % convert moment units from N.mm to N.m
%% construct the instrument penetration into the trocar
Ls_min = 0.160;
Ls_max = L - 0.01;
numofCycles = 5;
Ls_period = (size(ft,1))/numofCycles;
Ls = Ls_min+(Ls_max-Ls_min)*(sin((1:size(ft,1))'/2/Ls_period*2*pi)).^2;
%% calculate the sensor normalized output
shaftObj = shaftObj.ft_to_csF(ft,Ls,Lp);
shaftObj = shaftObj.F_to_Nn();
%% plot the forces
figure(1)
set(gcf, 'WindowState', 'maximized');
subplot(3,2,1:2)
plot(Ls,'b','linewidth',2)
plotsettings('sequence Number','Ls (m)')

subplot(3,2,3);
plot(ft(:,1),'b','linewidth',2)
hold on
plot(ft(:,2),'r','linewidth',2)
plotsettings('sequence Number','F(N)')
legend('Fx','Fy');

subplot(3,2,4);
plot(ft(:,4),'b','linewidth',2)
hold on
plot(ft(:,5),'r','linewidth',2)
plotsettings('sequence Number','M(N.m)')
legend('Mx','My');


subplot(3,2,5);
plot(ft(:,3),'b','linewidth',2)
plotsettings('sequence Number','Fz(N)')

subplot(3,2,6);
plot(ft(:,6),'b','linewidth',2)
plotsettings('sequence Number','Mz(N.m)')

figure(2)
plot(shaftObj.Nn,'linewidth',2);
set(gcf, 'WindowState', 'maximized');
plotsettings('sequence Number','N_{ni} (dmls)')
legend('i=1','i=2','i=3','i=4','i=5','i=6')
%% solve the optimization function with specified bounds
cx = 6*shaftObj.E*shaftObj.Ixx/shaftObj.ks;
cy = 6*shaftObj.E*shaftObj.Iyy/shaftObj.ks;

lb = [-inf*ones(36,1);0.6*L;0.2*cx;0.2*cy;1.5*Lp;-L];
ub = [ inf*ones(36,1);1.4*L; 10*cx; 10*cy;1.5*Lp;L];

randomLsOffset = rand(1)*20/1000;
encoderLs = Ls + randomLsOffset;

f = @(x) upNonLinModel(x,encoderLs,shaftObj.Nn,shaftObj.tpF);

options = optimoptions(@lsqnonlin,'MaxIterations',1000,'Algorithm',...
'trust-region-reflective','Display','iter','MaxFunctionEvaluations',inf);
[optimal_x,resnorm] = lsqnonlin(f,rand(41,1),lb,ub,options);

id_A = inv(reshape(optimal_x(1:36),6,6));
%% display the results
disp("----------------------------------------------------------------")
fprintf("%-20s%-20s%-20s\n",'Variable','True Value','Identified Value');
fprintf("%-20s%-20f%-20f\n",'L',L,optimal_x(37));
fprintf("%-20s%-20f%-20f\n",'cx',cx,optimal_x(38));
fprintf("%-20s%-20f%-20f\n",'cy',cy,optimal_x(39));
fprintf("%-20s%-20f%-20f\n",'Lp',Lp,optimal_x(40));
fprintf("%-20s%-20f%-20f\n",'Ls_offset',-randomLsOffset,optimal_x(41));
disp("----------------------------------------------------------------")
disp("True A - Sensor calibration Matrix {F = A*Nn}")
disp(shaftObj.A)
disp("Identified A - Sensor calibration Matrix {F = A*Nn}")
disp(id_A)
%%
function ati = readForceData()
ati = table2array(readtable('atiData.txt'));
ati = ati-mean(ati(1:500,:));
ati(:,1:2) = -ati(:,1:2);
ati(:,3)= -ati(:,3);
ati(:,6) = ati(:,6);
end

function plotsettings(xlabelstring,ylabelstring)
grid on
set(gca,'linewidth',2);
set(gca,'fontsize',12);
xlabel(xlabelstring,'fontsize',14);
ylabel(ylabelstring,'fontsize',14);
end