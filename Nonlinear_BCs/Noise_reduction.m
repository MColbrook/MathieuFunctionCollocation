% This code produces figures 7 and 8 of the paper

clear; close all

Kvec=10.^(-2:0.03:2);
L=length(Kvec);
h=@(x) 0.6*(0.2969*sqrt((x+1)/2)-0.126*((x+1)/2)-0.3516*((x+1)/2).^2+0.2843*((x+1)/2).^3-0.1015*((x+1)/2).^4);

PW_h1_lin_case1=0*Kvec;PW_h2_lin_case1=0*Kvec;
PW_h1_nonlin_case2=0*Kvec;PW_h2_nonlin_case2=0*Kvec;PW_h1_lin_case2=0*Kvec;PW_h2_lin_case2=0*Kvec;
PW_h1_nonlin_case3=0*Kvec;PW_h2_nonlin_case3=0*Kvec;PW_h1_lin_case3=0*Kvec;PW_h2_lin_case3=0*Kvec;
PW_h1_nonlin_case4=0*Kvec;PW_h2_nonlin_case4=0*Kvec;PW_h1_lin_case4=0*Kvec;PW_h2_lin_case4=0*Kvec;

% case 1
K=0;
for j=1:L
    PW_h1_lin_case1(j)=PW_comp(Kvec(j),@(x) 0*x+0.012,K,0,0);
    PW_h2_lin_case1(j)=PW_comp(Kvec(j),h,K,0,0);
end
%%
% case 2
K=2.7*10^(-9); beta=0.14;
for j=1:L
    PW_h1_nonlin_case2(j)=PW_comp(Kvec(j),@(x) 0*x+0.012,K,beta,1);
    PW_h1_lin_case2(j)=PW_comp(Kvec(j),@(x) 0*x+0.012,K,0,1);
    PW_h2_nonlin_case2(j)=PW_comp(Kvec(j),h,K,beta,1);
    PW_h2_lin_case2(j)=PW_comp(Kvec(j),h,K,0,1);
end

% case 3
K=5.72*10^(-11); beta=0.5;
for j=1:L
    PW_h1_nonlin_case3(j)=PW_comp(Kvec(j),@(x) 0*x+0.012,K,beta,0);
    PW_h1_lin_case3(j)=PW_comp(Kvec(j),@(x) 0*x+0.012,K,0,0);
    PW_h2_nonlin_case3(j)=PW_comp(Kvec(j),h,K,beta,0);
    PW_h2_lin_case3(j)=PW_comp(Kvec(j),h,K,0,0);
end

% case 4
K=3.65*10^(-12); beta=0.613;
for j=1:L
    PW_h1_nonlin_case4(j)=PW_comp(Kvec(j),@(x) 0*x+0.012,K,beta,0);
    PW_h1_lin_case4(j)=PW_comp(Kvec(j),@(x) 0*x+0.012,K,0,0);
    PW_h2_nonlin_case4(j)=PW_comp(Kvec(j),h,K,beta,0);
    PW_h2_lin_case4(j)=PW_comp(Kvec(j),h,K,0,0);
end

%%
h = plot(1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11);
c = get(h,'Color');
close all
%% h=0.012 plot
figure
semilogx(Kvec,PW_h1_lin_case1,'color',c{1},'linewidth',1.5)
hold on
semilogx(Kvec,PW_h1_lin_case2,'color',c{2},'linewidth',1.5)
semilogx(Kvec,PW_h1_nonlin_case2,'color',c{2},'linestyle','--','linewidth',1.5)
semilogx(Kvec,PW_h1_lin_case3,'color',c{3},'linewidth',1.5)
semilogx(Kvec,PW_h1_nonlin_case3,'color',c{3},'linestyle','--','linewidth',1.5)
semilogx(Kvec,PW_h1_lin_case4,'color',c{4},'linewidth',1.5)
semilogx(Kvec,PW_h1_nonlin_case4,'color',c{4},'linestyle','--','linewidth',1.5)
xlabel('$k_0$','interpreter','latex','fontsize',14)
ylabel('Scattered Far-field Sound','interpreter','latex','fontsize',14)
legend({'Case 1','Case 2, linear','Case 2, non-linear','Case 3, linear','Case 3, non-linear','Case 4, linear','Case 4, non-linear'},'interpreter','latex','fontsize',14,'location','northwest')
axis([0.01,100,-55,10])
%% h variable case
figure
semilogx(Kvec,PW_h2_lin_case1,'color',c{1},'linewidth',1.5)
hold on
semilogx(Kvec,PW_h2_lin_case2,'color',c{2},'linewidth',1.5)
semilogx(Kvec,PW_h2_nonlin_case2,'color',c{2},'linestyle','--','linewidth',1.5)
semilogx(Kvec,PW_h2_lin_case3,'color',c{3},'linewidth',1.5)
semilogx(Kvec,PW_h2_nonlin_case3,'color',c{3},'linestyle','--','linewidth',1.5)
semilogx(Kvec,PW_h2_lin_case4,'color',c{4},'linewidth',1.5)
semilogx(Kvec,PW_h2_nonlin_case4,'color',c{4},'linestyle','--','linewidth',1.5)
xlabel('$k_0$','interpreter','latex','fontsize',14)
ylabel('Scattered Far-field Sound','interpreter','latex','fontsize',14)
legend({'Case 1','Case 2, linear','Case 2, non-linear','Case 3, linear','Case 3, non-linear','Case 4, linear','Case 4, non-linear'},'interpreter','latex','fontsize',14,'location','northwest')
axis([0.01,100,-55,10])

%%
figure
semilogx(Kvec,PW_h1_lin_case1-PW_h1_lin_case2,'k','linewidth',2)
hold on
semilogx(Kvec,PW_h1_lin_case1-PW_h1_nonlin_case2,'--k','linewidth',2)
semilogx(Kvec,PW_h2_lin_case1-PW_h2_lin_case2,'k','linewidth',0.5)
semilogx(Kvec,PW_h2_lin_case1-PW_h2_nonlin_case2,'--k','linewidth',0.5)

XX=[4.1215 4.9458 6.5944 8.2430 9.8916 13.1888 16.4860 20.6075];
YY=[10 9.8000 8 6 3.5000 -0.4000 -4 -8.2000];

semilogx(XX(1:3),YY(1:3),'b.','markersize',20)
pp=semilogx(XX(4:end),YY(4:end),'.','markersize',20);

xlabel('$k_0$','interpreter','latex','fontsize',14)
ylabel('Noise Reduction','interpreter','latex','fontsize',14)
legend({'Linear, constant $h$','Non-linear, constant $h$','Linear, variable $h$','Non-linear, variable $h$','Experimental data'},'interpreter','latex','fontsize',14,'location','northwest')

axis([0.01,100,-10,40])
%%
figure
semilogx(Kvec,PW_h1_lin_case1-PW_h1_lin_case3,'k','linewidth',2)
hold on
semilogx(Kvec,PW_h1_lin_case1-PW_h1_nonlin_case3,'--k','linewidth',2)
semilogx(Kvec,PW_h2_lin_case1-PW_h2_lin_case3,'k','linewidth',0.5)
semilogx(Kvec,PW_h2_lin_case1-PW_h2_nonlin_case3,'--k','linewidth',0.5)

XX=[2.1981,2.7477,5.4953,8.2430,13.7383,16.4860,19.2337,21.9813,24.7290,27.4767,34.3458];
YY=[2.6543,0.7407,-4.6296,-4.2284,1.7284,4.1605,4.5679,5.1852,4.6605,5.0926,2.7778];

semilogx(XX([1,6:end]),YY([1,6:end]),'b.','markersize',20)
semilogx(XX(2:5),YY(2:5),'.','markersize',20)

xlabel('$k_0$','interpreter','latex','fontsize',14)
ylabel('Noise Reduction','interpreter','latex','fontsize',14)
legend({'Linear, constant $h$','Non-linear, constant $h$','Linear, variable $h$','Non-linear, variable $h$','Experimental data'},'interpreter','latex','fontsize',14,'location','northwest')
axis([0.01,100,-5,20])

%%
function PW = PW_comp(K0,h,K,beta,TT)
rho=1.225;
mu=1.81*10^(-5);
c0=343;
len=0.075; %plate semichord
Cf= beta; %intertial resistance, also called beta in various places

Grf=rho*Cf*c0*sqrt(K)/mu;
Grn=rho*K*c0/(mu*len);

P{1}.a0=-1; % left endpoint
P{1}.b0=1; % right endpoint
P{1}.M=min(1500,max(100,round(K0*3))); % number of M functions used in expansion
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

x=-1:1:1;
THETA=0:0.001:2*pi; 
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
    P{1}.C0 = @(x) (x>0.6)*K0^2;
end

[~,d_s,~,~,~] = Nonlinear_Porous_Plate(K0,phi_I,P,[],THETA,x);
PW=log10(sum(abs(d_s).^2)*0.001)*10;
end