clear
clc
% close all
%% set up the parameters

K0=25; % frequency

% incident wave
theta=pi/8;
phi_I.fun=@(x,y) exp(-1i*K0*(x*cos(theta)+y*sin(theta)))*0.1;
phi_I.dx=@(x,y) -1i*K0*cos(theta)*exp(-1i*K0*(x*cos(theta)+y*sin(theta)))*0.1; % x derivative
phi_I.dy=@(x,y) -1i*K0*sin(theta)*exp(-1i*K0*(x*cos(theta)+y*sin(theta)))*0.1; % y derivative

% grid of points where we compute the scattered field
Y=-0.99:0.02:1;
X=-2:0.02:2;
[X0,Y0] = meshgrid(X,Y); % WARNING: reduce size for speed
Z1=X0(:)+1i*Y0(:);
x=-1:0.001:1; % points on plate where we compute displacement etc. (rescaled to lie in lie in [-1,1])
THETA=0:0.01:2*pi; % angles where we compute the far-field directivity

% plate parameters
P=cell(1,1); % cell capturing parameters for each plate

P{1}.N=200; % number of Chebyshev polynomials used in expansion of eta. Set to zero if you want no elasticity (only porous plate)
P{1}.M=200; % number of M functions in expansion of jump in field across plate
P{1}.Mplot1=P{1}.M; % number of M functions for grid plot
P{1}.Mplot2=P{1}.M; % number of M functions for far field plot

P{1}.bc_type = 1; % 1=clamped-clamped,2=clamped-free,3=free-clamped,4=free-free
P{1}.a0=-1; % left end-point of plate
P{1}.b0=1; % left end-point of plate

% physical parameters in paper (x rescaled to lie in lie in [-d,d])
P{1}.B = @(x) 1000-200*sin(x);
P{1}.Bd = @(x) -200*cos(x); % first derivative of B
P{1}.Bdd = @(x) 200*sin(x); % second derivative of B
P{1}.R = @(x) 0*x+0.01;
P{1}.alphaH = @(x) 0.03+0*x;

nu = 0.25; % Poisson ratio
rho = 1.23; % a standard air density
c_0=343; % speed of sound
omega=K0*c_0; % angular frequency
m = 1; % mass per unit area

% define constants for each plate
for j=1:length(P)
    P{j}.B0 = @(x) -(1-P{j}.alphaH(x)).*m*omega^2;
    P{j}.B4 = @(x) (1-P{j}.alphaH(x)).*(1-2*P{j}.alphaH(x).*nu./(1-nu)).*P{j}.B(x);
    P{j}.B1 = @(x) 0*x;
    P{j}.B2 = @(x) (1-P{j}.alphaH(x)).*(1-2*P{j}.alphaH(x).*nu./(1-nu)).*P{j}.Bdd(x);
    P{j}.B3 = @(x) 2*(1-P{j}.alphaH(x)).*(1-2*P{j}.alphaH(x).*nu./(1-nu)).*P{j}.Bd(x);
    
    P{j}.nu = nu;
    P{j}.rho = rho;
    P{j}.c_0=c_0;
    P{j}.omega=omega;
    P{j}.m = m;
end



%% run the computation
[phi_s,d_s,phi_i,PHI_jump,ETA,ETA_d] = Plate_scatter(K0,phi_I,P,Z1,THETA,x); % see Plate_scatter.m for explanation of variables

%% plot the results
if isempty(Z1)==0
    figure
    Q=reshape(phi_s,length(Y),length(X));
    imagesc(X,Y,real(Q));
    set(gca,'YDir','normal');
    set(gca,'color','r');
    colormap('jet')
    colorbar
    axis tight
    title('Scattered Field')
    hold on
    for j=1:length(P)
        plot([real(P{j}.a0),real(P{j}.b0)],[imag(P{j}.a0),imag(P{j}.b0)],'k','linewidth',1.5)
    end
    axis equal
    
    Qi=reshape(phi_i,length(Y),length(X));
    figure
    imagesc(X,Y,real(Qi+Q));
    set(gca,'YDir','normal');
    set(gca,'color','r');
    colormap('jet')
    colorbar
    axis tight
    title('Total Field')
    hold on
    for j=1:length(P)
        plot([real(P{j}.a0),real(P{j}.b0)],[imag(P{j}.a0),imag(P{j}.b0)],'k','linewidth',1.5)
    end
    axis equal
end


figure
polarplot(THETA,abs(d_s))
title('Far-Field Directivity')

for j=1:length(P)
    figure
    plot(x,real(PHI_jump{j}))
    hold on
    plot(x,imag(PHI_jump{j}))
    title('Jump in Field across Plate')
    
    figure
    plot(x,real(ETA{j}))
    hold on
    plot(x,imag(ETA{j}))
    title('Plate Displacement')
    
    figure
    plot(x,real(ETA_d{j}))
    hold on
    plot(x,imag(ETA_d{j}))
    title('Derivative of Plate Displacement')
end
