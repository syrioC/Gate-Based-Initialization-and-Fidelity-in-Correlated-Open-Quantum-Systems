
sigx = [0 1;1 0];
sigp = [0 1;0 0];
sigm = [0 0;1 0];
sigz = [1 0;0 -1];
sigy = [0 -1i;1i 0];
alpha=0.01;
s=1;
kT=0;
xi=2;
pulse_length = 50;
theta = pi/2;
Delta = 1;
omegac =1;
[rho_taup2] = Finited_Gate_Withingate(alpha, s, kT, xi, pulse_length, theta, omegac);
dt = 0.001;
td = round(pulse_length / dt);
[rho_list1, rho_list2, rho_list3, t_evoN]  = Finite_Gate_Postgate(rho_taup2, alpha, s, kT, dt, theta, xi, Delta, td, nt, omegac);
sx1_list = zeros(1,numel(t_evoN));
sy1_list = zeros(1,numel(t_evoN));
sz1_list = zeros(1,numel(t_evoN));
sx2_list = zeros(1,numel(t_evoN));
sy2_list = zeros(1,numel(t_evoN));
sz2_list = zeros(1,numel(t_evoN));
sx3_list = zeros(1,numel(t_evoN));
sy3_list = zeros(1,numel(t_evoN));
sz3_list = zeros(1,numel(t_evoN));
for idx = 1:numel(t_evoN)
    rho1 = reshape(rho_list1(idx,:,:),2,2);
    rho2 = reshape(rho_list2(idx,:,:),2,2);
    rho3 = reshape(rho_list3(idx,:,:),2,2);
    sx1_list(1,idx) = real(trace(sigx*rho1))/2;
    sy1_list(1,idx) = real(trace(sigy*rho1))/2;
    sz1_list(1,idx) = real(trace(sigz*rho1))/2;
    sx2_list(1,idx) = real(trace(sigx*rho2))/2;
    sy2_list(1,idx) = real(trace(sigy*rho2))/2;
    sz2_list(1,idx) = real(trace(sigz*rho2))/2;
    sx3_list(1,idx) = real(trace(sigx*rho3))/2;
    sy3_list(1,idx) = real(trace(sigy*rho3))/2;
    sz3_list(1,idx) = real(trace(sigz*rho3))/2;
end

figure(1); 
plot(t_evoN,sx1_list); hold on
plot(t_evoN,sx2_list); hold on
plot(t_evoN,sx3_list);
title('Sx versus t')
legend('dynamically prepared', 'factorized initial condition', 'markovian')

figure(2); 
plot(t_evoN,sy1_list); hold on
plot(t_evoN,sy2_list); hold on
plot(t_evoN,sy3_list);
title('Sy versus t')
legend('dynamically prepared', 'factorized initial condition', 'markovian')
figure(3); 
plot(t_evoN,sz1_list); hold on
plot(t_evoN,sz2_list); hold on
plot(t_evoN,sz3_list);
title('Sz versus t')
legend('dynamically prepared', 'factorized initial condition', 'markovian')