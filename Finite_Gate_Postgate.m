function [rho_list, rho_list2, rho_list3, t_evoN]  = Finite_Gate_Postgate(rho0, alpha, s, kT, dt, theta, xi, Delta, td, nt, omegac)

if nargin < 11, omegac = 1; end
if nargin < 10, nt = 2^18; end
if nargin < 9, td = td; end
if nargin < 8, Delta = 1; end
sigx=[0 1;1 0];
sigp = [0 1;0 0];
sigm = [0 0;1 0];
sigz=[1 0;0 -1];
sigy=[0 -1i;1i 0];

%parameters
%rho_0 = [1 0; 0 0];

omegac=1;
nt = 2^18;


pulse_length = dt*td;


Delta = 1;
omegap = theta/(pulse_length);

Uc = expm(-1i*theta/2 *sigx);
rho_0 = rho0;

ID = [1 0; 0 1];
H_0 = - Delta*sigz/2;

[T,bcf]=BCFkT(nt*dt,nt,omegac,s,kT);
[Grp,~,Tr]=timedDkT(T,bcf,Delta);
[Grm,~,~]=timedDkT(T,bcf,-Delta);
[Grpp,~,~]=timedDkT(T,bcf,Delta+omegap);
[Grpm,~,~]=timedDkT(T,bcf,Delta-omegap);
[Grmp,~,~]=timedDkT(T,bcf,-Delta+omegap);
[Grmm,~,~]=timedDkT(T,bcf,-Delta-omegap);
[Gr0p,~,~]=timedDkT(T,bcf,omegap);
[Gr0m,~,~]=timedDkT(T,bcf,-omegap);
[Gr0,~,~]=timedDkT(T,bcf,0);

Ur=[1 0 0 1;
        0 1 1i 0;
        0 1 -1i 0;
        1 0 0 -1]/sqrt(2);

Dt0=-1i*(kron(eye(2),H_0)-kron(H_0.',eye(2)));
Dt0r=real(Ur'*Dt0*Ur);

DBR_static_IP = zeros(4, 4);


DBR_static_IP(2,2) = -2*xi^2*real(Gr0(end)) -1/2 * (real(Grp(end))+real(Grm(end)));
DBR_static_IP(2,3) =  1/2 * (imag(Grp(end))-imag(Grm(end)));
DBR_static_IP(3,2) =  1/2 * (imag(Grm(end))-imag(Grp(end)));
DBR_static_IP(3,3) = -2*xi^2*real(Gr0(end)) -1/2 * (real(Grp(end))+real(Grm(end)));
DBR_static_IP(4,1) =  -real(Grm(end))+real(Grp(end));
DBR_static_IP(4,4) =  -real(Grm(end))-real(Grp(end));

% Da = 1*alpha/2*(Ur'*reshape(DBR_static_IP(:,:),4,4)*Ur)+ Dt0;
% [Vr,~]=eigs(Da);
% rho_0 = reshape(Vr(:,4),2,2);
% rho_0 = rho_0/trace(rho_0)
% Uc = expm(-1i*theta/2 *sigx);
% rho_0 = Uc*rho_0*Uc'; 


A_0 = (sigx+xi*sigz)/2;
rho_t = rho_0;
rho_list = zeros(nt/4-td-2,2,2);
rho_list(1,:,:) = rho_t;
rho_list2 = zeros(nt/4-td-2,2,2);
rho_list2(1,:,:) = rho_t;
rho_list3 = zeros(nt/4-td-2,2,2);
rho_list3(1,:,:) = rho_t;
h = 4*dt;
i0=1;
t_evo = Tr(1:end-1);
Slambda1=zeros(numel(Tr),2,2);
Slambdainf=zeros(numel(Tr),2,2);
Slambda2 =zeros(numel(Tr),2,2);
Slambda3 =zeros(numel(Tr),2,2);
Slambda4 =zeros(numel(Tr),2,2);
Slambda =zeros(numel(Tr),2,2);
for i= 1:numel(Tr)
    Slambda1(i,:,:) = alpha*(sigp*Grp(i)/2 + sigm*Grm(i)/2 + xi*sigz*Gr0(i)/2); 
end
for i= 1:numel(Tr)
    Slambdainf(i,:,:) = alpha*(sigp*Grp(end)/2 + sigm*Grm(end)/2 + xi*sigz*Gr0(end)/2); 
end

for i= 1:numel(Tr)-td
    Phit = -theta*Tr(i)/pulse_length;
    part1 = 1/2 * (sigp*deltaGamma(i+td,i, Grp)+ sigm*deltaGamma(i+td,i, Grm));
    part2 = xi/4 * (sigp*(exp(1i*Phit)*deltaGamma(i+td,i, Grpp)-exp(-1i*Phit)*deltaGamma(i+td,i, Grpm))...
        +sigm*(exp(1i*Phit)*deltaGamma(i+td,i, Grmp)-exp(-1i*Phit)*deltaGamma(i+td,i, Grmm)));
    part3 = xi/4*sigz*(exp(1i*Phit)*deltaGamma(i+td,i, Gr0p)+exp(-1i*Phit)*deltaGamma(i+td,i, Gr0m));
    Slambda3(i,:,:) = alpha*(part1+part2+part3);
end

for i= 1:numel(Tr)-td
    Phi = -theta*Tr(i)/pulse_length;
    u = Tr(i+td);
    cz = xi/4*(exp(1j*Phi)*deltaGamma(i+td,i, Gr0p)+exp(-1j*Phi)*deltaGamma(i+td,i, Gr0m))...
        -1/8*(exp(1j*(Phi+Delta*u))*deltaGamma(i+td,i, Grpp)+exp(-1j*(Phi+Delta*u))*deltaGamma(i+td,i, Grmm))...
        +1/8*(exp(1j*(Phi-Delta*u))*deltaGamma(i+td,i, Grmp)+exp(-1j*(Phi-Delta*u))*deltaGamma(i+td,i, Grpm));
    cplus = 1/4*(exp(2*1j*u*Delta)*deltaGamma(i+td,i, Grm)+deltaGamma(i+td,i, Grp)...
        -exp(2*1j*u*Delta)/2*(exp(1j*Phi)*deltaGamma(i+td,i, Grmp)+exp(-1j*Phi)*deltaGamma(i+td,i, Grmm))...
        +1/2*(exp(1j*Phi)*deltaGamma(i+td,i, Grpp)+exp(-1j*Phi)*deltaGamma(i+td,i, Grpm))...
        +xi*exp(1j*Delta*u)*(exp(1j*Phi)*deltaGamma(i+td,i, Gr0p)-exp(-1j*Phi)*deltaGamma(i+td,i, Gr0m)));
    cminus = 1/4*(exp(-2*1j*u*Delta)*deltaGamma(i+td,i, Grp)+deltaGamma(i+td,i, Grm)...
        -exp(-2*1j*u*Delta)/2*(exp(1j*Phi)*deltaGamma(i+td,i, Grpp)+exp(-1j*Phi)*deltaGamma(i+td,i, Grpm))...
        +1/2*(exp(1j*Phi)*deltaGamma(i+td,i, Grmp)+exp(-1j*Phi)*deltaGamma(i+td,i, Grmm))...
        +xi*exp(-1j*Delta*u)*(exp(-1j*Phi)*deltaGamma(i+td,i, Gr0m)-exp(1j*Phi)*deltaGamma(i+td,i, Gr0p)));
    Slambda4(i,:,:) = alpha*(cz*sigz+cplus*sigp+cminus*sigm);
end

for i= 1:numel(Tr)
    Unit = expm(-1j*H_0*(Tr(i)-Tr(1)));
    UCX = Unit*Uc*Unit';
    Slambda2(i,:,:) = UCX*(reshape(Slambdainf(end,:,:),2,2)-reshape(Slambda1(i,:,:),2,2))*UCX';
end
for i= 1:numel(Tr)-td
    Slambda(i,:,:) = Slambda1(i,:,:)+Slambda2(i+td,:,:) +Slambda4(i,:,:) ;
end


for j=1:nt/4-td/2-2
    rho_list(j,:,:) = rho_t;
    
    % k1 = drhodt(i0,lambda,rho_t,reshape(At(i0,:,:),2,2));
    % k2= drhodt(i0+1,lambda,rho_t+h*k1/2,reshape(At(i0+1,:,:),2,2));
    % k3= drhodt(i0+1,lambda,rho_t+h*k2/2,reshape(At(i0+1,:,:),2,2));
    % k4 = drhodt(i0+2,lambda,rho_t+h*k3,reshape(At(i0+2,:,:),2,2));

    k1 = drhodtSanalytical(i0,Slambda,rho_t,A_0,H_0,t_evo);
    k2= drhodtSanalytical(i0+1,Slambda,rho_t+h*k1/2,A_0,H_0,t_evo);
    k3= drhodtSanalytical(i0+1,Slambda,rho_t+h*k2/2,A_0,H_0,t_evo);
    k4 = drhodtSanalytical(i0+2,Slambda,rho_t+h*k3,A_0,H_0,t_evo);

    i0=i0+2;
    rho_t=rho_t+h*(k1+2*k2+2*k3+k4)/6;
end
rho_t = rho_0;
i0=1;
for j=1:nt/4-td/2-2
    rho_list2(j,:,:) = rho_t;
    
    % k1 = drhodt(i0,lambda,rho_t,reshape(At(i0,:,:),2,2));
    % k2= drhodt(i0+1,lambda,rho_t+h*k1/2,reshape(At(i0+1,:,:),2,2));
    % k3= drhodt(i0+1,lambda,rho_t+h*k2/2,reshape(At(i0+1,:,:),2,2));
    % k4 = drhodt(i0+2,lambda,rho_t+h*k3,reshape(At(i0+2,:,:),2,2));

    k1 = drhodtSanalytical(i0,Slambda1,rho_t,A_0,H_0,t_evo);
    k2= drhodtSanalytical(i0+1,Slambda1,rho_t+h*k1/2,A_0,H_0,t_evo);
    k3= drhodtSanalytical(i0+1,Slambda1,rho_t+h*k2/2,A_0,H_0,t_evo);
    k4 = drhodtSanalytical(i0+2,Slambda1,rho_t+h*k3,A_0,H_0,t_evo);

    i0=i0+2;
    rho_t=rho_t+h*(k1+2*k2+2*k3+k4)/6;
end
rho_t = rho_0;
i0=1;
for j=1:nt/4-td/2-2
    rho_list3(j,:,:) = rho_t;
    
    % k1 = drhodt(i0,lambda,rho_t,reshape(At(i0,:,:),2,2));
    % k2= drhodt(i0+1,lambda,rho_t+h*k1/2,reshape(At(i0+1,:,:),2,2));
    % k3= drhodt(i0+1,lambda,rho_t+h*k2/2,reshape(At(i0+1,:,:),2,2));
    % k4 = drhodt(i0+2,lambda,rho_t+h*k3,reshape(At(i0+2,:,:),2,2));

    k1 = drhodtSanalytical(i0,Slambdainf,rho_t,A_0,H_0,t_evo);
    k2= drhodtSanalytical(i0+1,Slambdainf,rho_t+h*k1/2,A_0,H_0,t_evo);
    k3= drhodtSanalytical(i0+1,Slambdainf,rho_t+h*k2/2,A_0,H_0,t_evo);
    k4 = drhodtSanalytical(i0+2,Slambdainf,rho_t+h*k3,A_0,H_0,t_evo);

    i0=i0+2;
    rho_t=rho_t+h*(k1+2*k2+2*k3+k4)/6;
end
rhoN = numel(rho_list(:,1,1));
t_evoN = 2*t_evo(1:rhoN);
end

function deltaG = deltaGamma(upper_idx,lower_idx, Gr)
    deltaG = Gr(upper_idx)-Gr(lower_idx);
end

function [drho]=drhodtSanalytical(i0,Slambda,rho_t,A_0,H_0,t_evo)
    drho = drhodt(i0,Slambda,rho_t,A_0)-1i*((H_0)*rho_t-rho_t*(H_0));
end

