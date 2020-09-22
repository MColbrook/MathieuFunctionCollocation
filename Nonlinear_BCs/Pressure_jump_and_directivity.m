% This code produces figures 9 and 10 of the paper

clear; close all
h=@(x) 0.6*(0.2969*sqrt((x+1)/2)-0.126*((x+1)/2)-0.3516*((x+1)/2).^2+0.2843*((x+1)/2).^3-0.1015*((x+1)/2).^4);


% case 3
K=5.72*10^(-11); beta=0.5;

[d_s,PHI_jump]=PW_comp(0.1,@(x) 0*x+0.012,K,beta,0);
D_nonlin_h1_1=d_s; p_nonlin_h1_1=PHI_jump{1};
[d_s,PHI_jump]=PW_comp(0.1,@(x) 0*x+0.012,K,0,0);
D_lin_h1_1=d_s; p_lin_h1_1=PHI_jump{1};

[d_s,PHI_jump]=PW_comp(1,@(x) 0*x+0.012,K,beta,0);
D_nonlin_h1_2=d_s; p_nonlin_h1_2=PHI_jump{1};
[d_s,PHI_jump]=PW_comp(1,@(x) 0*x+0.012,K,0,0);
D_lin_h1_2=d_s; p_lin_h1_2=PHI_jump{1};

[d_s,PHI_jump]=PW_comp(10,@(x) 0*x+0.012,K,beta,0);
D_nonlin_h1_3=d_s; p_nonlin_h1_3=PHI_jump{1};
[d_s,PHI_jump]=PW_comp(10,@(x) 0*x+0.012,K,0,0);
D_lin_h1_3=d_s; p_lin_h1_3=PHI_jump{1};



[d_s,PHI_jump]=PW_comp(0.1,h,K,beta,0);
D_nonlin_h2_1=d_s; p_nonlin_h2_1=PHI_jump{1};
[d_s,PHI_jump]=PW_comp(0.1,h,K,0,0);
D_lin_h2_1=d_s; p_lin_h2_1=PHI_jump{1};

[d_s,PHI_jump]=PW_comp(1,h,K,beta,0);
D_nonlin_h2_2=d_s; p_nonlin_h2_2=PHI_jump{1};
[d_s,PHI_jump]=PW_comp(1,h,K,0,0);
D_lin_h2_2=d_s; p_lin_h2_2=PHI_jump{1};

[d_s,PHI_jump]=PW_comp(10,h,K,beta,0);
D_nonlin_h2_3=d_s; p_nonlin_h2_3=PHI_jump{1};
[d_s,PHI_jump]=PW_comp(10,h,K,0,0);
D_lin_h2_3=d_s; p_lin_h2_3=PHI_jump{1};

%%
figure
subplot(2,3,1)
plot(-1:0.01:1,real(p_nonlin_h1_1),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_nonlin_h1_1),'linewidth',2,'linestyle','--')
title('Non-linear, $k_0=0.1$','interpreter','latex','fontsize',10)

subplot(2,3,2)
plot(-1:0.01:1,real(p_nonlin_h1_2),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_nonlin_h1_2),'linewidth',2,'linestyle','--')
title('Non-linear, $k_0=1$','interpreter','latex','fontsize',10)

subplot(2,3,3)
plot(-1:0.01:1,real(p_nonlin_h1_3),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_nonlin_h1_3),'linewidth',2,'linestyle','--')
title('Non-linear, $k_0=10$','interpreter','latex','fontsize',10)

subplot(2,3,4)
plot(-1:0.01:1,real(p_lin_h1_1),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_lin_h1_1),'linewidth',2,'linestyle','--')
title('Linear, $k_0=0.1$','interpreter','latex','fontsize',10)

subplot(2,3,5)
plot(-1:0.01:1,real(p_lin_h1_2),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_lin_h1_2),'linewidth',2,'linestyle','--')
title('Linear, $k_0=1$','interpreter','latex','fontsize',10)

subplot(2,3,6)
plot(-1:0.01:1,real(p_lin_h1_3),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_lin_h1_3),'linewidth',2,'linestyle','--')
title('Linear, $k_0=10$','interpreter','latex','fontsize',10)


%%
figure
subplot(2,3,1)
plot(-1:0.01:1,real(p_nonlin_h2_1),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_nonlin_h2_1),'linewidth',2,'linestyle','--')
title('Non-linear, $k_0=0.1$','interpreter','latex','fontsize',10)

subplot(2,3,2)
plot(-1:0.01:1,real(p_nonlin_h2_2),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_nonlin_h2_2),'linewidth',2,'linestyle','--')
title('Non-linear, $k_0=1$','interpreter','latex','fontsize',10)

subplot(2,3,3)
plot(-1:0.01:1,real(p_nonlin_h2_3),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_nonlin_h2_3),'linewidth',2,'linestyle','--')
title('Non-linear, $k_0=10$','interpreter','latex','fontsize',10)

subplot(2,3,4)
plot(-1:0.01:1,real(p_lin_h2_1),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_lin_h2_1),'linewidth',2,'linestyle','--')
title('Linear, $k_0=0.1$','interpreter','latex','fontsize',10)

subplot(2,3,5)
plot(-1:0.01:1,real(p_lin_h2_2),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_lin_h2_2),'linewidth',2,'linestyle','--')
title('Linear, $k_0=1$','interpreter','latex','fontsize',10)

subplot(2,3,6)
plot(-1:0.01:1,real(p_lin_h2_3),'linewidth',2)
hold on
plot(-1:0.01:1,imag(p_lin_h2_3),'linewidth',2,'linestyle','--')
title('Linear, $k_0=10$','interpreter','latex','fontsize',10)


%%

figure

subplot(2,1,1)
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_nonlin_h1_1))*10,':k','linewidth',2);
hold on
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_nonlin_h1_2))*10,'--k','linewidth',1);
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_nonlin_h1_3))*10,'k','linewidth',2);
axis([0,1,-40,-10])
xlabel('$\theta/(2\pi)$','interpreter','latex','fontsize',14)
legend({'$k_0=0.1$','$k_0=1$','$k_0=10$'},'interpreter','latex','fontsize',10,'location','best')

subplot(2,1,2)
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_lin_h1_1))*10,':k','linewidth',2);
hold on
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_lin_h1_2))*10,'--k','linewidth',1);
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_lin_h1_3))*10,'k','linewidth',2);
axis([0,1,-40,-10])
xlabel('$\theta/(2\pi)$','interpreter','latex','fontsize',14)
legend({'$k_0=0.1$','$k_0=1$','$k_0=10$'},'interpreter','latex','fontsize',10,'location','best');


%%

figure

subplot(2,1,1)
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_nonlin_h2_1))*10,':k','linewidth',2);
hold on
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_nonlin_h2_2))*10,'--k','linewidth',1);
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_nonlin_h2_3))*10,'k','linewidth',2);
axis([0,1,-40,-10])
xlabel('$\theta/(2\pi)$','interpreter','latex','fontsize',14)
legend({'$k_0=0.1$','$k_0=1$','$k_0=10$'},'interpreter','latex','fontsize',10,'location','best')
subplot(2,1,2)
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_lin_h2_1))*10,':k','linewidth',2);
hold on
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_lin_h2_2))*10,'--k','linewidth',1);
plot((0:0.001:2*pi)/(2*pi),log10(abs(D_lin_h2_3))*10,'k','linewidth',2);
axis([0,1,-40,-10])
xlabel('$\theta/(2\pi)$','interpreter','latex','fontsize',14)
legend({'$k_0=0.1$','$k_0=1$','$k_0=10$'},'interpreter','latex','fontsize',10,'location','best')





%%
function [d_s,PHI_jump] = PW_comp(K0,h,K,beta,TT)
rho=1.225;
mu=1.81*10^(-5);
c0=343;
len=0.075; %plate semichord
Cf= beta; %intertial resistance, also called beta in various places

Grf=rho*Cf*c0*sqrt(K)/mu;
Grn=rho*K*c0/(mu*len);

P{1}.a0=-1; % left endpoint
P{1}.b0=1; % right endpoint
P{1}.M=min(1500,max(100,round(K0*3)));%min(1600,max(100,round(K0*4))); % number of M functions used in expansion
P{1}.N=P{1}.M; % number of M functions used in expansion
P{1}.Mplot1=P{1}.M;
P{1}.Mplot2=P{1}.M;

% physical variable parameters (x rescaled to lie in lie in [-d,d])
if K>0
    P{1}.C0 = @(x) K0^(2)+0*x;
    P{1}.C1 = @(x) 1i*K0*h(x)/Grn;
    P{1}.C2 = @(x) 1i*Grf*h(x)*K0^2/Grn;
else
    P{1}.C0 = @(x) 0*x;
    P{1}.C1 = @(x) 1+0*x;
    P{1}.C2 = @(x) 0*x;
end
    
x=-1:0.01:1;
THETA=0:0.001:2*pi; % need small spacing of grid for large k0

%quadrupole
x00=0.99;
y00=0.1;


Z0 = x00 + 1i.*y00;
RR = @(x,y) sqrt((x-x00).^2+(y-y00).^2);
QUAD1 = @(x,y) 0.0144*(1i*K0.^2./(4.*RR(x,y).^2)).*( besselh(2,K0*RR(x,y)).*(x-real(Z0)).*(y-imag(Z0)));


phi_I.fun = @(x,y) QUAD1(x,y);
phi_I.dx = @(x,y) (QUAD1(x+0.00000001,y)-QUAD1(x-0.00000001,y))/0.00000002;
phi_I.dy = @(x,y) (QUAD1(x,y+0.00000001)-QUAD1(x,y-0.00000001))/0.00000002;

P{1}.SCALE=abs(P{1}.C1(0)); % the algorithm computes SCALE*Eta_A then rescales
if TT==1
    P{1}.C0 = @(x) (x>0.6)*K0^2;%(tanh((x-0.8)*10)+1)/2.*K0^2;
end

[~,d_s,~,PHI_jump,~] = Nonlinear_Porous_Plate(K0,phi_I,P,[],THETA,x);
end