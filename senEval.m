clc;clear all;close all
%%
H = 50/1000;
HprimeList = (35)/1000;
dz = 21/1000;
rs = 14.5/1000;

E = 44.31E9;
ri = 6.8/2/1000;
ro = 8.43/2/1000;
Ixx = 1/4*pi*(ro^4-ri^4);
Iyy = 1/4*pi*(ro^4-ri^4);
Jzz = 1/2*pi*(ro^4-ri^4);
Area = pi*(ro^2-ri^2);
G = 10E9;

Sw = 1/1000;
Cg = 0.1/1000;
Cp = (Sw-Cg)/2;

for i=1:numel(HprimeList)
    Hprime = HprimeList(i);
    Lambda = zeros(6);
    
    Lambda(1,1) = H^3/3/E/Iyy+Hprime*H^2/2/E/Iyy;
    Lambda(1,5) = H^2/2/E/Iyy;
    
    Lambda(2,2) = H^3/3/E/Ixx+Hprime*H^2/2/E/Ixx;
    Lambda(2,4) = -H^2/2/E/Ixx;
    
    Lambda(3,3) = H/Area/E;
    
    Lambda(4,2) = -H^2/2/E/Ixx - Hprime*H/E/Ixx;
    Lambda(4,4) = H/E/Ixx;
    
    Lambda(5,1) = H^2/2/E/Iyy + Hprime*H/E/Iyy;
    Lambda(5,5) = H/E/Iyy;
    
    Lambda(6,6) = H/G/Jzz;
    
    dxx = 0.0/1000;
    dyy = 0.0/1000;
    dzz = 0.0/1000;
    
    Psi = zeros(6);
    Psi(1,2) = 1;
    Psi(1,4) = dz;
    Psi(1,6) = rs+dxx;
    
    Psi(2,3) = 1;
    Psi(2,4) = sqrt(3)/2*rs+dyy;
    Psi(2,5) = -1/2*rs-dxx;
    
    Psi(3,1) = -sqrt(3)/2;
    Psi(3,2) = -1/2;
    Psi(3,4) = 1/2*(dzz-dz);
    Psi(3,5) = sqrt(3)/2*(dz-dzz);
    Psi(3,6) = rs-dxx/2+dyy*sqrt(3)/2;
    
    Psi(4,3) = 1;
    Psi(4,4) = dyy;
    Psi(4,5) = rs-dxx;
    
    Psi(5,1) = sqrt(3)/2;
    Psi(5,2) = -1/2;
    Psi(5,4) = 1/2*(dzz-dz);
    Psi(5,5) = sqrt(3)/2*(dzz-dz);
    Psi(5,6) = rs-1/2*dxx-sqrt(3)/2*dyy;
    
    Psi(6,3) = 1;
    Psi(6,4) = -sqrt(3)/2*rs+dyy;
    Psi(6,5) = -rs/2-dxx;
    
    theta = 5/180*pi;
      
    C = 2*Psi*Lambda/Cp*diag([1/10,1/10,1/10,1/0.4,1/0.4,1/0.4]);
   
    condnum = cond(C,inf)
    [U,S,V] = svd(C)
    x(:,i) = diag(S);
end
%%
counter = 1;
figure;
set(gcf,'position',[84, 558, 1752, 420]);
subplot(1,4,1)
plot(HprimeList*1000,x(1,:),'b','linewidth',2)
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 2;
grid on
xlabel('H^''(mm)','fontsize',16,'fontweight','b');
ylim([0.98*min(x(1,:)) 1.02*max(x(1,:))]);
ylabel('\Sigma[1,1]','fontsize',16,'fontweight','b');

subplot(1,4,2)
plot(HprimeList*1000,x(2,:),'b','linewidth',2)
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 2;
grid on
xlabel('H^''(mm)','fontsize',16,'fontweight','b');
ylim([0.98*min(x(2,:)) 1.02*max(x(2,:))]);
ylabel('\Sigma[2,2]','fontsize',16,'fontweight','b');

subplot(1,4,3)
plot(HprimeList*1000,x(4,:),'b','linewidth',2)
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 2;
grid on
xlabel('H^''(mm)','fontsize',16,'fontweight','b');
ylim([0.98*min(x(4,:)) 1.02*max(x(4,:))]);
ylabel('\Sigma[4,4]','fontsize',16,'fontweight','b');

subplot(1,4,4)
plot(HprimeList*1000,x(6,:),'b','linewidth',2)
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 2;
grid on
xlabel('H^''(mm)','fontsize',16,'fontweight','b');
ylim([0.98*min(x(6,:)) 1.02*max(x(6,:))]);
ylabel('\Sigma[6,6]','fontsize',16,'fontweight','b');
print('condNumstudy','-djpeg','-r600')