function [rhofinal] = Finited_Gate_Withingate(alpha, s, kT, xi, tp, theta, omegac)
    td = 4096; 
    nt = 8*td;
    % Pauli matrices
    sigx = [0 1;1 0];
    sigp = [0 1;0 0];
    sigm = [0 0;1 0];
    sigz = [1 0;0 -1];
    sigy = [0 -1i;1i 0];
    ID = eye(2);
    pulse_length = tp;
    % parameters
    omegap = theta/pulse_length;
    
    Delta = 1;
    dt = pulse_length/td;

    [T,bcf] = BCFkT(nt*dt,nt,omegac,s,kT);
    [Grp,~,Tr] = timedDkT(T,bcf,Delta);
    [Grm,~,~] = timedDkT(T,bcf,-Delta);
    [Grpp,~,~] = timedDkT(T,bcf,Delta+omegap);
    [Grpm,~,~] = timedDkT(T,bcf,Delta-omegap);
    [Grmp,~,~] = timedDkT(T,bcf,-Delta+omegap);
    [Grmm,~,~] = timedDkT(T,bcf,-Delta-omegap);
    [Gr0p,~,~] = timedDkT(T,bcf,omegap);
    [Gr0m,~,~] = timedDkT(T,bcf,-omegap);
    [Gr0,~,~]  = timedDkT(T,bcf,0);

    Ur=[1 0 0 1;
        0 1 1i 0;
        0 1 -1i 0;
        1 0 0 -1]/sqrt(2);
    H_0 = -Delta*sigz/2;
    Dt0=-1i*(kron(eye(2),H_0)-kron(H_0.',eye(2)));
    Dt0r=real(Ur'*Dt0*Ur);
    
    DBR_static_SP = zeros(4, 4);
    
    DBR_static_SP(2,1) = -xi*(real(Grp(end))-real(Grm(end))) ;
    DBR_static_SP(2,2) = -2*xi^2*real(Gr0(end)) ;
    DBR_static_SP(2,3) =  0;
    DBR_static_SP(2,4) = xi*(real(Grp(end))+real(Grm(end))) ;
    DBR_static_SP(3,1) =  xi * (imag(Grm(end))+imag(Grp(end))-2*imag(Gr0(end)));
    DBR_static_SP(3,2) =   imag(Grm(end))-imag(Grp(end));
    DBR_static_SP(3,3) = -2*xi^2*real(Gr0(end)) - (real(Grp(end))+real(Grm(end)));
    DBR_static_SP(3,4) =  xi * (imag(Grm(end))-imag(Grp(end)));
    DBR_static_SP(4,1) =  -real(Grm(end))+real(Grp(end));
    DBR_static_SP(4,2) =  2*xi*real(Gr0(end));
    DBR_static_SP(4,4) =  -real(Grm(end))-real(Grp(end));
    
   Da = alpha/2*(reshape(DBR_static_SP(:,:),4,4))+ Dt0r;
   [Vr,Dr]=eigs(Da);
   rho_0 = Vr(1,4)*[1 0; 0 1]/2 + Vr(2,4)*sigx/2+Vr(3,4)*sigy/2+Vr(4,4)*[1 0; 0 -1]/2;
   %rho_0 = (rho_0 + rho_0') / 2;
   rho_0 = rho_0/trace(rho_0);


    H_0 = - Delta * sigz / 2;
    A_0 = (sigx + xi * sigz) / 2;

    

    h = 4 * dt;
    t_evo = Tr(1:end-1);
    i0 = 1;

    Slambda1 = zeros(numel(Tr),2,2);
    Slambdainf = zeros(numel(Tr),2,2);
    Slambda3 = zeros(numel(Tr),2,2);
    Slambda2 = zeros(numel(Tr),2,2);
    Slambda = zeros(numel(Tr),2,2);

    for i = 1:numel(Tr)
        Slambda1(i,:,:) = alpha*(sigp*Grp(i)/2 + sigm*Grm(i)/2 + xi*sigz*Gr0(i)/2); 
        Slambdainf(i,:,:) = alpha*(sigp*Grp(end)/2 + sigm*Grm(end)/2 + xi*sigz*Gr0(end)/2); 
    end

    for it = 1:numel(Tr)-1
        partz = xi/4*(Gr0p(it)+Gr0m(it)) - 1/8 * ...
            (exp(1i*Delta*t_evo(it))*(Grmp(it)-Grmm(it)) - exp(-1i*Delta*t_evo(it))*(Grpp(it)-Grpm(it)));
        partp = 1/4*(Grp(it) + 0.5*(Grpp(it)+Grpm(it))) ...
            + 0.25 * exp(2i*Delta*t_evo(it)) * (Grm(it) - 0.5*(Grmp(it)+Grmm(it))) ...
            + 0.25 * exp(1i*Delta*t_evo(it)) * (Gr0p(it)-Gr0m(it));
        partm = 1/4*(Grm(it) + 0.5*(Grmp(it)+Grmm(it))) ...
            + 0.25 * exp(-2i*Delta*t_evo(it)) * (Grp(it) - 0.5*(Grpp(it)+Grpm(it))) ...
            - 0.25 * exp(-1i*Delta*t_evo(it)) * (Gr0p(it)-Gr0m(it));
        Slambda3(it,:,:) = alpha * (partz*sigz + partp*sigp + partm*sigm);
    end

    for i = 1:numel(Tr)
        Unit = expm(-1i*H_0*(Tr(i)-Tr(1)));
        Uip = expm(-1i*omegap*sigx/2*(Tr(i)-Tr(1)));
        UCX = Unit * Uip * Unit';
        Slambda2(i,:,:) = UCX * (reshape(Slambdainf(end,:,:),2,2) - reshape(Slambda1(i,:,:),2,2)) * UCX';
    end

    for i = 1:numel(Tr)-td
        Slambda(i,:,:) = Slambda3(i,:,:) + Slambda2(i,:,:);
        %Slambda(i,:,:) = Slambdainf(i,:,:);
    end

    rho_t = rho_0;
    rho_list = zeros(nt/4 - td - 2, 2, 2);
    rhoref_list = zeros(nt/4 - td - 2, 2, 2);
    
    rho_list(1,:,:) = rho_t;
    rhoref_list(1,:,:) = rho_t;
    i0 = 1;
    for j = 1:nt/4 - td/2 - 2
        rho_list(j,:,:) = rho_t;
        k1 = drhodtSanalytical(i0,Slambda,rho_t,A_0,H_0,t_evo,omegap);
        k2 = drhodtSanalytical(i0+1,Slambda,rho_t + h*k1/2,A_0,H_0,t_evo,omegap);
        k3 = drhodtSanalytical(i0+1,Slambda,rho_t + h*k2/2,A_0,H_0,t_evo,omegap);
        k4 = drhodtSanalytical(i0+2,Slambda,rho_t + h*k3,A_0,H_0,t_evo,omegap);
        i0 = i0 + 2;
        rho_t = rho_t + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    end
    
    rhoN = numel(rho_list(:,1,1));
    t_evoN = 2*t_evo(1:rhoN);
    endidx = round(pulse_length/(t_evoN(2)-t_evoN(1)))+1;
    rho_tf = reshape(rho_list(endidx,:,:),2,2);
    rhofinal = rho_tf;
end

function [drho]=drhodtSanalytical(i0,Slambda,rho_t,A_0,H_0, t_evo, omegap)
    sigx=[0 1;1 0];
    Uint = expm(-1j*H_0*(t_evo(i0)-t_evo(1)));
    H_c_IP = omegap/2*sigx;
    H_c_SP = Uint*H_c_IP*Uint';
    drho = drhodt(i0,Slambda,rho_t,A_0)-1i*((H_c_SP+H_0)*rho_t-rho_t*(H_c_SP+H_0));
end


