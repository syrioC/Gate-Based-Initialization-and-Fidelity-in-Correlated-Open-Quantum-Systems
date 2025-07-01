
function [rho_opt, rho_base] = simulate_fidelity(lambda2, xi)
nt = 20000;
tp = 200;
a_opt = [-1.059484345450 0.439617195071 -0.118840413008 0.110510810142 -0.243876642586 -0.303546271785 0.376484113915];

tlist = linspace(0.01, tp, nt);
Delta = 1;
theta= pi/2;
sigx = [0 1; 1 0];  % sigma_x
sigz = [1 0; 0 -1]; 
sigy=[0 -1i;1i 0];
H_0 = -Delta*sigz/2;
A = (1*sigx+xi*sigz)/2;
Mlist = compute_M_list(a_opt, tlist, tp, Delta, A);
bMlist = compute_M_list(a_opt0, tlist, tp, Delta, A);
rho_0 = [1 0; 0 0];
rho_t = rho_0;
rho_list = zeros(nt/2,2,2);
brho_list = zeros(nt/2,2,2);
refrho_list = zeros(nt/2,2,2);
h = 2*(tlist(2)-tlist(1));
omegap = theta/tp;
i0=1;
rho_t = rho_0;
for j=1:nt/2-2
    omegap = df_s(tlist(i0),a_opt,tp);
    rho_list(j,:,:) = rho_t;
    k1 = drhodtSanalytical(i0,Mlist,rho_t,A,H_0,tlist,omegap);
    k2= drhodtSanalytical(i0+1,Mlist,rho_t+h*k1/2,A,H_0,tlist,omegap);
    k3= drhodtSanalytical(i0+1,Mlist,rho_t+h*k2/2,A,H_0,tlist,omegap);
    k4 = drhodtSanalytical(i0+2,Mlist,rho_t+h*k3,A,H_0,tlist,omegap);
    i0=i0+2;
    rho_t=rho_t+h*(k1+2*k2+2*k3+k4)/6;
end
i0=1;
rho_t = rho_0;
for j=1:nt/2-2
    omegap = df_s(tlist(i0),a_opt0,tp);
    brho_list(j,:,:) = rho_t;
    k1 = drhodtSanalytical(i0,bMlist,rho_t,A,H_0,tlist,omegap);
    k2= drhodtSanalytical(i0+1,bMlist,rho_t+h*k1/2,A,H_0,tlist,omegap);
    k3= drhodtSanalytical(i0+1,bMlist,rho_t+h*k2/2,A,H_0,tlist,omegap);
    k4 = drhodtSanalytical(i0+2,bMlist,rho_t+h*k3,A,H_0,tlist,omegap);
    i0=i0+2;
    rho_t=rho_t+h*(k1+2*k2+2*k3+k4)/6;
end

i0=1;
rho_t = rho_0;
for j=1:nt/2-2
    omegap = df_s(tlist(i0),a_opt0,tp);
    refrho_list(j,:,:) = rho_t;
    k1 = drhodtSanalytical(i0,0*bMlist,rho_t,A,H_0,tlist,omegap);
    k2= drhodtSanalytical(i0+1,0*bMlist,rho_t+h*k1/2,A,H_0,tlist,omegap);
    k3= drhodtSanalytical(i0+1,0*bMlist,rho_t+h*k2/2,A,H_0,tlist,omegap);
    k4 = drhodtSanalytical(i0+2,0*bMlist,rho_t+h*k3,A,H_0,tlist,omegap);
    i0=i0+2;
    rho_t=rho_t+h*(k1+2*k2+2*k3+k4)/6;
end
rhoN = numel(rho_list(:,1,1))-2;
t_evoN = 2*tlist(1:rhoN);
sxlist = zeros(1,rhoN);
sylist = zeros(1,rhoN);
szlist = zeros(1,rhoN);
bsxlist = zeros(1,rhoN);
bsylist = zeros(1,rhoN);
bszlist = zeros(1,rhoN);
refsxlist = zeros(1,rhoN);
refsylist = zeros(1,rhoN);
refszlist = zeros(1,rhoN);
rho_opt = reshape(rho_list(end-2,:,:),2,2);
rho_base = reshape(brho_list(end-2,:,:),2,2);
fidlist = zeros(1,rhoN);
bfidlist = zeros(1,rhoN);
Uc = expm(-1i*theta/2 *sigx);
ref_rho = Uc*rho_0*Uc'; 
for j=1:rhoN 
    rho = reshape(rho_list(j,:,:),2,2);
    brho = reshape(brho_list(j,:,:),2,2);
    Uint = expm(-1i*t_evoN(j) *H_0);
    refrho = Uint*ref_rho*Uint';
    sxlist(j) = real(trace(rho*sigx/2));
    sylist(j) = real(trace(rho*sigy/2));
    szlist(j) = real(trace(rho*sigz/2));
    fidlist(j) = 1-0.5*sum(svd(rho-refrho));
    bsxlist(j) = real(trace(brho*sigx/2));
    bsylist(j) = real(trace(brho*sigy/2));
    bszlist(j) = real(trace(brho*sigz/2));
    bfidlist(j) = 1-0.5*sum(svd(brho-refrho));
    %fidlist(j) = 1-0.5*sum(svd(rho-rho2));
end
fid_opt = fidlist(end);
fid_base = bfidlist(end);
end 
function U = Ux(s, t, Delta, fval)
    sigx = [0 1; 1 0];
    sigz = [1 0; 0 -1];
    H0 = -Delta/2 * sigz;
    U = expm(-1i * H0 * t) * expm(-1i * fval/2 * sigx) * expm(1i * H0 * t);
end
function R = Rmatrix(s, t, Delta, A, fval)
    sigz = [1 0; 0 -1];
    H0 = -Delta/2 * sigz;
    U = Ux(s, t, Delta, fval);
    U_H = expm(-1i * H0 * s);
    R = U * U_H * A * U_H' * U';
end
function c = Corr(s)
    lambda2 = 1e-5;
    omegac=1;
    c = 2*lambda2*omegac^2*2/(1+1i*omegac*s)^(1+1/2);
end

function val = f_s(s, a, tp)
    N = length(a);
    val = (pi/2) * s / tp;
    for n = 1:N
        val = val + a(n) * sin(n * pi * s / tp);
    end
end

function val = df_s(s, a, tp)
    N = length(a);
    val = (pi/2)/ tp;
    for n = 1:N
        val = val + n * pi / tp *a(n) * cos(n * pi * s / tp);
    end
end

function Mlist = compute_M_list(a_opt, tlist, tp, Delta, A)
    Ntime = length(tlist);
    Mlist = zeros(2, 2, Ntime);
    sigz = [1 0; 0 -1];
    smax_tail = 400;
    
    for k = 1:Ntime
        t = tlist(k);
        Ns = 500;
        s_vals = linspace(0, t, Ns);
        Rsum = zeros(2,2);

        % --- First integral: from s = 0 to t ---
        for j = 1:Ns
            s = s_vals(j);
            fval = f_s(s, a_opt, tp);
            R = Rmatrix(s, t, Delta, A, fval);
            Rsum = Rsum + Corr(s) * R;
        end
        M1 = Rsum * (t / Ns);

        % --- Second integral: from s = t to smax_tail ---
        Nt_tail = 500;
        s_tail = linspace(t, smax_tail, Nt_tail);
        tail_sum = zeros(2,2);
        H0 = -Delta/2 * sigz;
        for j = 1:Nt_tail
            s = s_tail(j);
            U = expm(-1i * H0 * s);
            tail_sum = tail_sum + Corr(s) * (U * A * U');
        end
        tail_avg = tail_sum * ((smax_tail - t)/Nt_tail);

        % --- Apply Ux(t) on the tail ---
        fval_t = f_s(t, a_opt, tp);
        Uxt = Ux(t, tp, Delta, fval_t);
        M2 = Uxt * tail_avg * Uxt';

        % --- Final M(t) ---
        Mlist(:, :, k) = M1 + M2;
    end
end


function [drho]=drhodtSanalytical(i0,Slambda,rho_t,A_0,H_0, t_evo, omegap)
    sigx=[0 1;1 0];
    Uint = expm(-1j*H_0*(t_evo(i0)-t_evo(1)));
    H_c_IP = omegap/2*sigx;
    H_c_SP = Uint*H_c_IP*Uint';
    drho = drhodt(i0,Slambda,rho_t,A_0)-1i*((H_c_SP+H_0)*rho_t-rho_t*(H_c_SP+H_0));
end



