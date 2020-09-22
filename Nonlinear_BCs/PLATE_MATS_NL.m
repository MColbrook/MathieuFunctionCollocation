function [MAT1,MAT2_A,MAT2_B,MAT2_C,M0,M0_dx,M0_dy,M0_far,M0_jump,eta] = PLATE_MATS_NL(K0,P,Z1,Z2,THETA,x)
% Constructs the submatrices needed for collocation on multiple plates 
%%% INPUTS %%%
% P:     struct of plate variables
% Z1:    points where we compute the scattered field
% Z2:    points where we compute scattered field derivatives for collocation with other plates
% THETA: angles where we compute far-field
% x:     points along plate where we compute jump, eta etc. rescaled to lie in [-1,1]

%%% OUTPUTS %%%
% MAT1:      collocation matrix for kinematic condition on plate
% MAT2_A:    collocation matrix A for nonlinear BCs
% MAT2_B:    collocation matrix B for nonlinear BCs
% MAT2_C:    collocation matrix C for nonlinear BCs
% M0:        value matrix at points Z1
% M0_dx:     value matrix of x derivative at points Z2
% M0_dy:     value matrix of y derivative at points Z2
% M0_far:    value matrix of far-field at angles THETA
% M0_jump:   value matrix of jump in field at x
% eta:       value matrix of eta_a

% Author: Matthew Colbrook, http://www.damtp.cam.ac.uk/user/mjc249/home.html

%% Define various parameters in terms of P
M=P.M;
N=P.N;
Mplot1=P.Mplot1; % number of Mathieu functions for evaluation at Z1
Mplot2=P.Mplot2; % number of Mathieu functions for evaluation at THETA

a0=P.a0; % left end-point in complex form
b0=P.b0; % right end-point in complex form
ANG=angle(b0-a0); % angle of plate relative to usual coordinates
d=abs(a0-b0)/2; % scaling parameter for Sturm-Liouville problem
Q=d^2*K0^2/4; % rescaled Sturm-Liouville parameter

%% Pre-compute spectral info for Mathieu functions
[~,B0] = Eig_problem(0,Q,round((M+1)/2)); % these are used to construct Mathieu functions (much quicker to do once and save)
[~,B1] = Eig_problem(1,Q,round((M+1)/2));

B=B0(:,1:round((M+1)/2));
X0=zeros(1,round((M+1)/2));

for j=1:round((M+1)/2)
    X0=X0+(-1)^j*mtimes((  -(sqrt(Q)).*(besselasym(j-2,j+1,0,sqrt(Q))-besselasym(j,j+1,0,sqrt(Q)))/2+(sqrt(Q)).*(besselasym(j-1,j,0,sqrt(Q))-besselasym(j-1,j+2,0,sqrt(Q)))/2 ... 
                       +(sqrt(Q)).*(besselasym(j,j-1,0,sqrt(Q))-besselasym(j+2,j-1,0,sqrt(Q)))/2-(sqrt(Q)).*(besselasym(j+1,j-2,0,sqrt(Q))-besselasym(j+1,j,0,sqrt(Q)))/2),B(j,:));
end

B=B1(:,1:round((M+1)/2));
X1=zeros(1,round((M+1)/2));

for j=1:round((M+1)/2)
    X1=X1-(-1)^j*mtimes((  -(sqrt(Q)).*(besselasym(j-2,j,0,sqrt(Q))-besselasym(j,j,0,sqrt(Q)))/2+(sqrt(Q)).*(besselasym(j-1,j-1,0,sqrt(Q))-besselasym(j-1,j+1,0,sqrt(Q)))/2 ...
                      +(sqrt(Q)).*(besselasym(j-1,j-1,0,sqrt(Q))-besselasym(j+1,j-1,0,sqrt(Q)))/2-(sqrt(Q)).*(besselasym(j,j-2,0,sqrt(Q))-besselasym(j,j,0,sqrt(Q)))/2   ),B(j,:));
end

save('Eig_dat','B0','B1','X0','X1')

%% Set up collocation matrix 1 (kinematic condition)
xpts1=d*cos(pi*(2*(1:M)-1)/(2*M));

M1=M_per(1,round((M+1)/2),ceil(M/2),acos(xpts1(:)/d));
M2=M_per(0,round((M+1)/2),ceil(M/2),acos(xpts1(:)/d));

MAT1=[M1,M2];
MAT1=MAT1(:,1:M);

T=zeros(N,length(xpts1));
T(1,:)=1;
T(2,:)=xpts1/d;
for j=3:N
    T(j,:)=2*xpts1/d.*T(j-1,:)-T(j-2,:);
end
T1=T;
for j=1:length(xpts1)
    T1(:,j)=-P.C0(xpts1(j))*sqrt(d^2-xpts1(j)^2)*T(:,j);
end
MAT1=[MAT1,transpose(T1)/P.SCALE];

%% Set up collocation matrix 2
xpts2=d*cos(pi*(2*(1:N)-1)/(2*N));

M1=M_per(1,round((M+1)/2),ceil(M/2),acos(xpts2(:)/d));
M2=M_per(0,round((M+1)/2),ceil(M/2),acos(xpts2(:)/d));
MH1=M_hankel(1,Q,round((M+1)/2),ceil(M/2),0);
MH2=M_hankel(0,Q,round((M+1)/2),ceil(M/2),0);

MAT2_A=2*[M1,M2].*[MH1,MH2];
MAT2_A=MAT2_A(:,1:M);

T=zeros(N,length(xpts2));
T(1,:)=1;
T(2,:)=xpts2/d;
for j=3:N
    T(j,:)=2*xpts2/d.*T(j-1,:)-T(j-2,:);
end

T1=T;T2=T;
for j=1:N
    T1(j,:)=-P.C1(xpts2).*T1(j,:);
    T2(j,:)=-P.C2(xpts2).*T2(j,:);
end


MAT2_B=[0*MAT2_A,transpose(T2)]/P.SCALE;
MAT2_C=[0*MAT2_A,transpose(T)]/P.SCALE;
MAT2_A=[MAT2_A,transpose(T1)/P.SCALE];
clear T1 T2

%% Compute eta_a
T=zeros(N,length(x));
T(1,:)=1;
T(2,:)=x;
for j=3:N
    T(j,:)=2*x.*T(j-1,:)-T(j-2,:);
end
eta=transpose(T);

%% Compute scattered field
if isempty(Z1)==0
    Z1=(Z1-(a0+b0)/2)*exp(-1i*ANG);
    Zhat=acosh(Z1(:)/d);
    TAU=imag(Zhat);
    NU=real(Zhat);

    M1=M_per(1,round((M+1)/2),ceil(Mplot1/2),TAU);
    M2=M_per(0,round((M+1)/2),ceil(Mplot1/2),TAU);

    MH1=M_hankel(1,Q,round((M+1)/2),ceil(Mplot1/2),NU);
    MH2=M_hankel(0,Q,round((M+1)/2),ceil(Mplot1/2),NU);

    M0=[M1,M2].*[MH1,MH2];
    M0=M0(:,1:Mplot1);
else
    M0=[];
end


if isempty(Z2)==0
    Z2=(Z2-(a0+b0)/2)*exp(-1i*ANG);
    Zhat=acosh(Z2(:)/d);
    TAU=imag(Zhat);
    NU=real(Zhat);

    M1=M_per(1,round((M+1)/2),ceil(M/2),TAU);
    M2=M_per(0,round((M+1)/2),ceil(M/2),TAU);
    MH1=M_hankel(1,Q,round((M+1)/2),ceil(M/2),NU);
    MH2=M_hankel(0,Q,round((M+1)/2),ceil(M/2),NU);

    M1_d=M_per_d(1,round((M+1)/2),ceil(M/2),TAU);
    M2_d=M_per_d(0,round((M+1)/2),ceil(M/2),TAU);
    MH1_d=M_hankel_d(1,Q,round((M+1)/2),ceil(M/2),NU);
    MH2_d=M_hankel_d(0,Q,round((M+1)/2),ceil(M/2),NU);


    M0_dtau=[M1_d,M2_d].*[MH1,MH2];M0_dtau=M0_dtau(:,1:M);
    M0_dnu=[M1,M2].*[MH1_d,MH2_d];M0_dnu=M0_dnu(:,1:M);

    % convert to rotated x and y
    M0_dxp=(repmat(sinh(NU).*cos(TAU)./((sinh(NU).*cos(TAU)).^2+(cosh(NU).*sin(TAU)).^2),1,M)/d).*M0_dnu...
        -(repmat(cosh(NU).*sin(TAU)./((sinh(NU).*cos(TAU)).^2+(cosh(NU).*sin(TAU)).^2),1,M)/d).*M0_dtau;

    M0_dyp=(repmat(cosh(NU).*sin(TAU)./((sinh(NU).*cos(TAU)).^2+(cosh(NU).*sin(TAU)).^2),1,M)/d).*M0_dnu...
        +(repmat(sinh(NU).*cos(TAU)./((sinh(NU).*cos(TAU)).^2+(cosh(NU).*sin(TAU)).^2),1,M)/d).*M0_dtau;

    % convert to unrotated coordinates
    M0_dx=cos(ANG)*M0_dxp-sin(ANG)*M0_dyp;
    M0_dy=sin(ANG)*M0_dxp+cos(ANG)*M0_dyp;

else
    M0_dx=[];
    M0_dy=[];
end

%% Compute jump in field across plate
Z=x;
Zhat=acosh(Z(:));
TAU=imag(Zhat);
M1=M_per(1,round((M+1)/2),ceil(M/2),TAU);
M2=M_per(0,round((M+1)/2),ceil(M/2),TAU);
MH1=M_hankel(1,Q,round((M+1)/2),ceil(M/2),0);
MH2=M_hankel(0,Q,round((M+1)/2),ceil(M/2),0);

M0_jump=2*[M1,M2].*repmat([MH1,MH2],length(x),1);
M0_jump=M0_jump(:,1:M);

%% Compute far-field

if isempty(Z2)
    %% WARNING - this only works for one plate
    THETA=THETA-ANG;
    M1=M_per(1,round((M+1)/2),ceil(Mplot2/2),THETA(:));
    M2=M_per(0,round((M+1)/2),ceil(Mplot2/2),THETA(:));
    MH1=M_hankel_far(1,Q,ceil(Mplot2/2));
    MH2=M_hankel_far(0,Q,ceil(Mplot2/2));

    M0_far=[M1,M2].*repmat([MH1,MH2],length(THETA),1);
    M0_far=M0_far(:,1:Mplot2)*sqrt(d);
else
    Zfar=10000*exp(1i*THETA(:));
    Zfar=(Zfar-(a0+b0)/2)*exp(-1i*ANG);
    Zhat=acosh(Zfar(:)/d);
    TAU=imag(Zhat);
    NU=real(Zhat);

    M1=M_per(1,round((M+1)/2),ceil(Mplot2/2),TAU);
    M2=M_per(0,round((M+1)/2),ceil(Mplot2/2),TAU);

    MH1=M_hankel(1,Q,round((M+1)/2),ceil(Mplot2/2),NU);
    MH2=M_hankel(0,Q,round((M+1)/2),ceil(Mplot2/2),NU);

    M0_far=[M1,M2].*[MH1,MH2];
    M0_far=M0_far(:,1:Mplot2)*sqrt(10000);
end
    
    
end



%%
%%%%%%%% FUNCTIONS USED IN COMPUTATION %%%%%%%%

function [b,B] = Eig_problem(parity,q,N)
if parity==0
    A=diag(zeros(1,N-1)+q,1)+diag(zeros(1,N-1)+q,-1)+diag(4*(1:N).^2);
else
    A=diag(zeros(1,N-1)+q,1)+diag(zeros(1,N-1)+q,-1)+diag((2*(1:N)-1).^2);
    A(1,1)=A(1,1)-q;
end
[B,D] = eig((A));
b=diag(D);
end

function [X] = M_per(parity,N,m,tau)
load('Eig_dat','B0','B1')
if parity==0
    X=mtimes(sin(2*mtimes(tau(:),1:N)),B0(1:N,1:m));
else
    X=mtimes(sin(mtimes(tau(:),2*(1:N)-1)),B1(1:N,1:m));
end
end

function [X] = M_per_d(parity,N,m,tau)
load('Eig_dat','B0','B1')
if parity==0
    X=mtimes(2*mtimes(tau(:)*0+1,1:N).*cos(2*mtimes(tau(:),1:N)),B0(1:N,1:m));
else
    X=mtimes(mtimes(tau(:)*0+1,2*(1:N)-1).*cos(mtimes(tau(:),2*(1:N)-1)),B1(1:N,1:m));
end
end

function [X] = M_hankel(parity,q,N,m,nu)
X0=0;X1=0;
load('Eig_dat','B0','B1','X0','X1')
X=zeros(length(nu),m);
X0=X0(:,1:m);
X1=X1(:,1:m);
if parity==0
    for j=1:N
        X=X+(-1)^j*mtimes((besselasym(j-1,j+1,nu,sqrt(q))- besselasym(j+1,j-1,nu,sqrt(q))),B0(j,1:m));
    end
    X=X.*repmat(1./X0,length(nu),1);
else
    for j=1:N
        X=X-(-1)^j*mtimes((besselasym(j-1,j,nu,sqrt(q))- besselasym(j,j-1,nu,sqrt(q))),B1(j,1:m));
    end
    X=X.*repmat(1./X1,length(nu),1);
end

end

function [X] = M_hankel_d(parity,q,N,m,nu)
X0=0;X1=0;
load('Eig_dat','B0','B1','X0','X1')
X=zeros(length(nu),m);
X0=X0(:,1:m);
X1=X1(:,1:m);
if parity==0
    for j=1:N
        X=X+(-1)^j*mtimes((besselasym_d(j-1,j+1,nu,sqrt(q))- besselasym_d(j+1,j-1,nu,sqrt(q))),B0(j,1:m));
    end
    X=X.*repmat(1./X0,length(nu),1);
else
    for j=1:N
        X=X-(-1)^j*mtimes((besselasym_d(j-1,j,nu,sqrt(q))- besselasym_d(j,j-1,nu,sqrt(q))),B1(j,1:m));
    end
    X=X.*repmat(1./X1,length(nu),1);
end

end

function [X] = M_hankel_far(parity,q,m)
X0=0;X1=0;
load('Eig_dat','B0','B1','X0','X1')
X0=X0(:,1:m);
X1=X1(:,1:m);
if parity==0
    X=sqrt(2)./sqrt(pi*sqrt(q)*2).*exp(-1i*pi/4)*B0(1,1:m);
    X=X./X0;
else
    X=-sqrt(2)./sqrt(pi*sqrt(q)*2).*exp(1i*pi/4)*B1(1,1:m);
    X=X./X1;
end
end

function [X] = besselasym(ll,mm,nu,x)
% x is sqrt(Q)
M=5;
if x<0.1
    if max(ll,mm)<15
        X=besselj(ll,exp(-nu)*x).*besselh(mm,1,exp(nu)*x);
    else
        y2=0*nu;
        y3=0*nu;
        for j=0:M
            X=1+0*ll;
            for k=0:j
                X=X./(ll+k);
            end
            y2=y2+(-1)^j./(gamma((j+1))).*X.*(exp(-nu)*(x)/2).^(2*j);

            X=1+0*ll;
            if (mm-ll)-j<0
                for k=1:(j-(mm-ll))
                    X=X./(ll-k);
                end
            elseif (mm-ll)-j>0
                for k=0:((mm-ll)-j-1)
                    X=X.*(ll+k);
                end
            end
            y3=y3-1i/pi.*(exp(nu)*(x)/2).^(2*j)*((x)/2).^(-(mm-ll))/gamma(j+1).*X;
        end
        X=y2.*y3.*exp(-nu*(ll+mm));
    end
else
    if max(ll,mm)<max(105,round(1.5*x))
        X=besselj(ll,exp(-nu)*x).*besselh(mm,1,exp(nu)*x);
    else
        y2=0*nu;
        y3=0*nu;
        for j=0:M
            X=1+0*ll;
            for k=0:j
                X=X./(ll+k);
            end
            y2=y2+(-1)^j./(gamma((j+1))).*X.*(exp(-nu)*(x)/2).^(2*j);

            X=1+0*ll;
            if (mm-ll)-j<0
                for k=1:(j-(mm-ll))
                    X=X./(ll-k);
                end
            elseif (mm-ll)-j>0
                for k=0:((mm-ll)-j-1)
                    X=X.*(ll+k);
                end
            end
            y3=y3-1i/pi.*(exp(nu)*(x)/2).^(2*j)*((x)/2).^(-(mm-ll))/gamma(j+1).*X;
        end
        X=y2.*y3.*exp(-nu*(ll+mm));
    end
end
end

function [X] = besselasym_d(ll,mm,nu,x)

X=x*0.5*exp(-nu).*(besselasym(ll-1,mm,nu,x)-besselasym(ll+1,mm,nu,x))...
    +x*0.5*exp(nu).*(besselasym(ll,mm-1,nu,x)-besselasym(ll,mm+1,nu,x));

end


