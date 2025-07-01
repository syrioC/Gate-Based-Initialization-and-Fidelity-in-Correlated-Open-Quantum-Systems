sigx=[0 1;1 0];
sigp = [0 1;0 0];
sigm = [0 0;1 0];
sigz=[1 0;0 -1];
sigy=[0 -1i;1i 0];

%parameters

alpha = 0.02;
s= 1;
kT = 0.00/1.5;
nt = 2^18;
dt = 0.005;
td = 2;
nt= nt;
pulse_length = dt*td;
omegac = 1;
xi = 1;
theta =1*pi/2;
Delta = 1;

% rho_0 = [1 0; 0 0];
% Uc = expm(-1i*theta/2 *sigx);
% rho_0 = Uc*rho_0*Uc'; 

ID = [1 0; 0 1];
H_0 = - Delta*sigz/2;

[T,bcf]=BCFkT(nt*dt,nt,omegac,s,kT);
[Grp,~,Tr]=timedDkT(T,bcf,Delta);
[Grm,~,~]=timedDkT(T,bcf,-Delta);
[Gr0,~,~]=timedDkT(T,bcf,0);

dGrp(:) = Grp(end)-Grp(:);
dGrm(:) = Grm(end)-Grm(:);
dGr0(:) = Gr0(end)-Gr0(:);

figure(1);
plot(real(dGr0)); hold on;
plot(real(dGrp)); hold on;
plot(real(dGrm));
N = numel(Tr);

Ur=[1 0 0 1;
        0 1 1i 0;
        0 1 -1i 0;
        1 0 0 -1]/sqrt(2);

Dt0=-1i*(kron(eye(2),H_0)-kron(H_0.',eye(2)));
Dt0r=real(Ur'*Dt0*Ur);

DBR_static_IP = zeros(N, 4, 4);
DBR_rotate_IP = zeros(N, 4, 4);
for it = 1:N
    DBR_static_IP(it,2,2) = -2*xi^2*real(Gr0(it)) -1/2 * (real(Grp(it))+real(Grm(it)));
    DBR_static_IP(it,2,3) =  1/2 * (imag(Grp(it))-imag(Grm(it)));
    DBR_static_IP(it,3,2) =  1/2 * (imag(Grm(it))-imag(Grp(it)));
    DBR_static_IP(it,3,3) = -2*xi^2*real(Gr0(it)) -1/2 * (real(Grp(it))+real(Grm(it)));
    DBR_static_IP(it,4,1) =  -real(Grm(it))+real(Grp(it));
    DBR_static_IP(it,4,4) =  -real(Grm(it))-real(Grp(it));
end

for it = 1:N
    DBR_rotate_IP(it,2,1) = 2*xi^2*imag(dGr0(it)) + 1/2 * (imag(dGrp(it))+imag(dGrm(it)));
    DBR_rotate_IP(it,2,3) =  1/2 * (imag(dGrp(it))+imag(dGrm(it)));
    DBR_rotate_IP(it,3,1) =  1/2 * (real(dGrm(it))-real(dGrp(it)));
    DBR_rotate_IP(it,3,3) = -1/2 * (real(dGrm(it))+real(dGrp(it)));
    DBR_rotate_IP(it,3,4) = -2*xi^2*real(dGr0(it));
    DBR_rotate_IP(it,4,1) =  -1/2*(real(dGrm(it))-real(dGrp(it)));
    DBR_rotate_IP(it,4,4) =  -1/2*(real(dGrm(it))+real(dGrp(it)));
end

Da = 1*alpha/2*(Ur'*reshape(DBR_static_IP(end,:,:),4,4)*Ur)+ Dt0;
[Vr,R]=eigs(Da);
rho_0 = reshape(Vr(:,4),2,2);
rho_0 = rho_0/trace(rho_0)
% rho_0 = [1 0; 0 0];
Uc = expm(-1i*theta/2 *sigx);
rho_0 = Uc*rho_0*Uc'; 
rhon_0 = 0.5*[1; real(trace(sigx*rho_0));real(trace(sigy*rho_0));real(trace(sigz*rho_0))];
h = 4*dt;

Dt = alpha/2*(DBR_static_IP+DBR_rotate_IP);
Dtfa = alpha/2*(DBR_static_IP);
Dtma = alpha/2*(reshape(DBR_static_IP(end,:,:),4,4));


rhon_t = rhon_0;
rhon_list = zeros(nt/4-2,4);
rhon_list(1,:,:) = rhon_t;
rhonfa_list = zeros(nt/4-2,4);
rhonfa_list(1,:,:) = rhon_t;
rhonma_list = zeros(nt/4-2,4);
rhonma_list(1,:,:) = rhon_t;
i0=1;
for j=1:N/2-2
    rhon_list(j,:,:) = rhon_t;    
    k1 = reshape(Dt(i0,:,:),4,4)*rhon_t;
    k2= reshape(Dt(i0+1,:,:),4,4)*(rhon_t+h*k1/2);
    k3= reshape(Dt(i0+1,:,:),4,4)*(rhon_t+h*k2/2);
    k4 = reshape(Dt(i0+2,:,:),4,4)*(rhon_t+h*k3);
    i0=i0+2;
    rhon_t=rhon_t+h*(k1+2*k2+2*k3+k4)/6;
end

rhon_t = rhon_0;
i0=1;
for j=1:N/2-2
    rhonfa_list(j,:,:) = rhon_t;    
    k1 = reshape(Dtfa(i0,:,:),4,4)*rhon_t;
    k2= reshape(Dtfa(i0+1,:,:),4,4)*(rhon_t+h*k1/2);
    k3= reshape(Dtfa(i0+1,:,:),4,4)*(rhon_t+h*k2/2);
    k4 = reshape(Dtfa(i0+2,:,:),4,4)*(rhon_t+h*k3);
    i0=i0+2;
    rhon_t=rhon_t+h*(k1+2*k2+2*k3+k4)/6;
end

rhon_t = rhon_0;
i0=1;
for j=1:N/2-2
    rhonma_list(j,:,:) = rhon_t;    
    k1 = Dtma*rhon_t;
    k2= Dtma*(rhon_t+h*k1/2);
    k3= Dtma*(rhon_t+h*k2/2);
    k4 = Dtma*(rhon_t+h*k3);
    i0=i0+2;
    rhon_t=rhon_t+h*(k1+2*k2+2*k3+k4)/6;
end

t_evo = Tr(1:end-1);
rhoN = numel(rhon_list(:,1,1));
t_evoN = 2*t_evo(1:rhoN);
sxl = real(rhon_list(:,2));
syl = real(rhon_list(:,3));
szl = real(rhon_list(:,4));
sxlfa = real(rhonfa_list(:,2));
sylfa = real(rhonfa_list(:,3));
szlfa = real(rhonfa_list(:,4));
sxlma = real(rhonma_list(:,2));
sylma = real(rhonma_list(:,3));
szlma = real(rhonma_list(:,4));
figure(2);

plot(t_evoN,sxl); hold on
plot(t_evoN,syl); hold on
plot(t_evoN,szl); 
legend('sx','sy','sz')
xlim([0 200])

figure(3);
plot(t_evoN,sxl.^2+syl.^2); hold on
legend('coherence')
xlim([0 200])


IP_rholist = zeros(nt/4-2,2,2);
IP_rholistfa = zeros(nt/4-2,2,2);
IP_rholistma = zeros(nt/4-2,2,2);
sxlist = [];
sxfalist = [];
sxmalist = [];
sylist = [];
szlist = [];

for it = 1: nt/4-2
    IP_rholist(it,:,:) = 0.5*ID + sxl(it)*sigx + syl(it)*sigy+szl(it)*sigz;
    IP_rholistfa(it,:,:) = 0.5*ID + sxlfa(it)*sigx + sylfa(it)*sigy+szlfa(it)*sigz;
    IP_rholistma(it,:,:) = 0.5*ID + sxlma(it)*sigx + sylma(it)*sigy+szlma(it)*sigz;
end 


for it =1:nt/4-2
    iprho = reshape(IP_rholist(it,:,:),2,2);
    Uint = expm(-1i*H_0*(t_evoN(it)-t_evoN(1)));

    sprho = Uint*iprho*Uint';
    iprho = reshape(IP_rholistfa(it,:,:),2,2);
    sprhofa = Uint*iprho*Uint';
    iprho = reshape(IP_rholistma(it,:,:),2,2);
    sprhoma = Uint*iprho*Uint';

    sxlist(it) = real(trace(sprho*sigx/2));
    sxfalist(it) = real(trace(sprhofa*sigx/2));
    sxmalist(it) = real(trace(sprhoma*sigx/2));
    sylist(it) = real(trace(sprho*sigy/2));
    szlist(it) = real(trace(sprho*sigz/2));
end

figure(4);
plot(t_evoN,sxlist); hold on
plot(t_evoN,sylist); hold on
plot(t_evoN,szlist);

legend('sx','sy','sz')
xlim([0 600])


figure(5);
plot(t_evoN,sxlist); hold on
plot(t_evoN,sxfalist); hold on
plot(t_evoN,sxmalist);

legend('dpsx','fasx','masx')
xlim([0 600])