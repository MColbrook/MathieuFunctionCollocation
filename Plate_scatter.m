function [phi_s,d_s,phi_i,PHI_jump,ETA,ETA_d] = Plate_scatter(K0,phi_I,P,Z1,THETA,x)
% Performs the collocation method
%%% INPUTS %%%
% K0:    wavenumber
% phi_I: incident field and its derivatives
% P:     structure of plate variables
% Z1:    points where we compute scattered (and incident) field
% THETA: angles where we compute far-field
% x:     points along plate where we compute jump lie in [-1,1]

%%% OUTPUTS %%%
% phi_s:    scattered field at points Z1
% d_s:      far-field directivity at angles THETA
% phi_i:    incident field at points Z1
% PHI_jump: jumps in field across plates
% ETA:      eta on plates
% ETA_d:    eta_deriv on plates

% Author: Matthew Colbrook, http://www.damtp.cam.ac.uk/user/mjc249/home.html

%% Make sure everything is the correct shape
Z1=Z1(:);
THETA=THETA(:);
for j=1:length(P)
    if isfield(P{j}, 'N')==0
        P{j}.N=0;
    elseif (P{j}.N<5)&&(P{j}.N>0)
        P{j}.N=0;
        warning('If plate is elastic, number of Chebyshev polynomials for deformation expansion must be greater than 4 to include boundary conditions. Setting to empty expansion.')
    end
end

%% Compute collocation points for each plate
Z2=cell(length(P),1);
Z0=Z2;
for j=1:length(P)
    Z2{j}=transpose((P{j}.b0-P{j}.a0)/2*cos(pi*(2*(1:P{j}.M)-1)/(2*P{j}.M))+(P{j}.a0+P{j}.b0)/2); % collocation points in complex plane for kinematic
    Z0{j}=transpose(abs(P{j}.b0-P{j}.a0)/2*sqrt(1-(cos(pi*(2*(1:P{j}.M)-1)/(2*P{j}.M)).^2))); % factor when constructing the full matrix
end

%% Compute collocation matrices
MAT1=cell(length(P),1);MAT2=MAT1;M0=MAT1;M0_dx=MAT1;M0_dy=MAT1;M0_far=MAT1;M0_jump=MAT1;eta=MAT1;eta_d=MAT1;
for j=1:length(P)
    Zcol=[];
    for jj=1:length(P)
        if abs(j-jj)>0
            Zcol=[Zcol(:);Z2{jj}(:)];
        end
    end
    [MAT1{j},MAT2{j},M0{j},M0_dx{j},M0_dy{j},M0_far{j},M0_jump{j},eta{j},eta_d{j}] = PLATE_MATS(K0,P{j},Z1,Zcol,THETA,x);
end

Mtot=zeros(length(P)+1,1);
Ntot=zeros(length(P)+1,1);
for j=1:length(P)
    Mtot(j+1)=Mtot(j)+P{j}.M;
    Ntot(j+1)=Ntot(j)+P{j}.N;
end

%% Set up the full collocation matrix for the kinematic condition
A_kin=zeros(Mtot(end),Mtot(end)+Ntot(end));
RHS_kin=zeros(Mtot(end),1);
for j=1:length(P)
    % collocate at the P{j}.M points Z2{j}
    B=zeros(P{j}.M,Mtot(end)+Ntot(end));
    ANG=angle(P{j}.b0-P{j}.a0);
    B(:,(Mtot(j)+Ntot(j)+1):(Mtot(j+1)+Ntot(j+1)))=MAT1{j};
    for jj=1:length(P) % contribution of jjth plate
        if abs(jj-j)>0
            Mstart=0;
            for jjj=1:j-1
                if abs(jjj-jj)>0
                	Mstart=Mstart+P{jjj}.M;
                end
            end
            B(:,(Mtot(jj)+Ntot(jj)+1):(Mtot(jj+1)+Ntot(jj)))=...
                (-sin(ANG)*M0_dx{jj}(Mstart+1:Mstart+P{j}.M,:)+cos(ANG)*M0_dy{jj}(Mstart+1:Mstart+P{j}.M,:)).*repmat(Z0{j},1,P{jj}.M);
        end
    end
    A_kin(Mtot(j)+1:Mtot(j+1),:)=B;
    RHS_kin(Mtot(j)+1:Mtot(j+1))=-Z0{j}.*(-sin(ANG)*phi_I.dx(real(Z2{j}),imag(Z2{j}))+cos(ANG)*phi_I.dy(real(Z2{j}),imag(Z2{j})));
end

%% Set up the full collocation matrix for the beam equation
A_beam=zeros(Ntot(end),Mtot(end)+Ntot(end));
RHS_beam=zeros(Ntot(end),1);
for j=1:length(P)
    B=zeros(P{j}.N,Mtot(end)+Ntot(end));
    B(:,(Mtot(j)+Ntot(j)+1):(Mtot(j+1)+Ntot(j+1)))=MAT2{j};
    A_beam(Ntot(j)+1:Ntot(j+1),:)=B;
end

%% Set up full system and solve
A=[A_kin;A_beam];
RHS=[RHS_kin;RHS_beam];

sc = 1./sum(abs(A),2);
A = bsxfun(@times,sc,A);
RHS = bsxfun(@times,sc,RHS);
sc = 1./sum(abs(A),1);
A = bsxfun(@times,sc,A);
c=A\RHS;
c=transpose(bsxfun(@times,sc,transpose(c)));

Mcoeff=cell(length(P),1);
Ncoeff=cell(length(P),1);
for j=1:length(P)
    Mcoeff{j}=c((Mtot(j)+Ntot(j)+1):(Mtot(j+1)+Ntot(j)));
    Ncoeff{j}=c((Mtot(j+1)+Ntot(j)+1):(Mtot(j+1)+Ntot(j+1)));
end

%% compute the scattered fields
PHI_S=cell(length(P),1);
phi_s=0*Z1;
if isempty(Z1)==0
    for j=1:length(P)
        I=[1:ceil(P{j}.Mplot1/2),(ceil(P{j}.M/2)+1):(ceil(P{j}.M/2)+ceil(P{j}.Mplot1/2))];
        I=I(1:P{j}.Mplot1);
        PHI_S{j}=M0{j}*Mcoeff{j}(I);
        phi_s=phi_s+PHI_S{j};
    end
end

%% compute the incident field
phi_i=phi_I.fun(real(Z1(:)),imag(Z1(:)));

%% compute the directivity
D_S=cell(length(P),1);
d_s=0*THETA;
for j=1:length(P)
    I=[1:ceil(P{j}.Mplot2/2),(ceil(P{j}.M/2)+1):(ceil(P{j}.M/2)+ceil(P{j}.Mplot2/2))];
    I=I(1:P{j}.Mplot2);
    D_S{j}=M0_far{j}*Mcoeff{j}(I);
    d_s=d_s+D_S{j};
end

%% compute jumps across plates
PHI_jump=cell(length(P),1);
for j=1:length(P)
    PHI_jump{j}=M0_jump{j}*Mcoeff{j};
end

%% compute plate displacements and derivatives
ETA=cell(length(P),1);
ETA_d=cell(length(P),1);
for j=1:length(P)
    if P{j}.N>0
        ETA{j}=eta{j}*Ncoeff{j};
        ETA_d{j}=eta_d{j}*Ncoeff{j};
    else
        ETA{j}=0*x;
        ETA_d{j}=0*x;
    end
end

end