
function optimize_f_s2()
    tp = 100;
    Delta = 1;
    xi = 2;
    Nbasis = 7;
    Nt = 500;
    sigz = [1 0; 0 -1];
    A = sigz / 2;
    if isempty(gcp('nocreate'))
        parpool('local', 8);
    end
    a0 = 0.2 * randn(Nbasis, 1);
    opts = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'MaxIterations', 300, ...
        'UseParallel', true, ...
        'Algorithm', 'interior-point');
    [a_opt, fval] = fmincon( ...
        @(a) loss_total(a, tp, Delta, A, Nbasis, Nt), ...
        a0, [], [], [], [], [], [], ...
        @(a) nonlcon_endpoint_deriv(a, tp, Nbasis), ...
        opts);

    % ==== baseline loss ====
    a_baseline = zeros(Nbasis, 1);
    loss_baseline = loss_total(a_baseline, tp, Delta, A, Nbasis, Nt);
    reduction_pct = 100 * (loss_baseline - fval) / loss_baseline;

    s_vals = linspace(0, tp, 200);
    f_opt_vals = arrayfun(@(s) f_s(s, a_opt, tp), s_vals);
    f_base_vals = (pi/2) * s_vals / tp;

    figure;
    plot(s_vals, f_opt_vals, 'b-', 'LineWidth', 2); hold on;
    plot(s_vals, f_base_vals, 'r--', 'LineWidth', 2);
    xlabel('s'); ylabel('f(s)');
    title(sprintf('Optimized f(s) vs. Baseline  (Reduction = %.2f%%)', reduction_pct));
    legend('Optimized', 'Baseline', 'Location', 'SouthEast');
    ylim([0, pi/2]);
    grid on;

    fprintf('\nOptimized coefficients:\n');
    fprintf('%.12f  ', a_opt);
    fprintf('\n');

    % save('a_opt_result.txt', 'a_opt', '-ascii');
end




function cval = Corr(s)
    lambda2 = 1;
    omegac = 1;
    cval = 2 * lambda2 * omegac^2 * 2 / (1 + 1i * omegac * s)^(1.5);
end

function [c, ceq] = nonlcon_endpoint_deriv(a, tp, Nbasis)
    baseline = pi / (2 * tp);
    ceq = [df_s(0, a, tp) - 0;
           df_s(tp, a, tp) - 0];
    c = [];
end

function U = Ux(s, tp, Delta, fval)
    H0 = -Delta/2 * [1 0; 0 -1];
    U1 = expm(-1i * H0 * tp);
    U2 = expm(-1i * fval/2 * [0 1; 1 0]);
    U3 = expm(1i * H0 * tp);
    U = U1 * U2 * U3;
end

function R = Rmatrix(s, t, Delta, A, fval)
    H0 = -Delta/2 * [1 0; 0 -1];
    U = Ux(s, t, Delta, fval);
    U_H = expm(-1i * H0 * s);
    R = U * U_H * A * U_H' * U';
end

function val = f_s(s, a, tp)
    N = length(a);
    base = (pi/2) * s / tp;
    correction = 0;
    for n = 1:N
        correction = correction + a(n) * sin(n * pi * s / tp);
    end
    val = base + correction;
end

function val = df_s(s, a, tp)
    N = length(a);
    val = (pi/2) / tp;
    for n = 1:N
        val = val + a(n) * (n*pi/tp) * cos(n * pi * s / tp);
    end
end

function M = compute_M(t, a, tp, Delta, A)
    sigz = [1 0; 0 -1];
    Ns = 500;
    Nt_tail = 500;
    s_vals = linspace(0, t, Ns);
    s_tail = linspace(t, 400, Nt_tail);

    Rsum = zeros(2,2);
    for j = 1:Ns
        s = s_vals(j);
        fval = f_s(s, a, tp);
        R = Rmatrix(s, t, Delta, A, fval);
        Rsum = Rsum + Corr(s) * R;
    end
    term1 = Rsum * (t / Ns);

    % Tail part
    tail_sum = zeros(2,2);
    for j = 1:Nt_tail
        s = s_tail(j);
        U = expm(1i * (Delta/2) * sigz * s);
        tail_sum = tail_sum + Corr(s) * (U * A * U');
    end
    tail_avg = tail_sum * ((max(s_tail) - t) / Nt_tail);

    % Apply Ux
    fval_t = f_s(t, a, tp);
    Uxt = Ux(t, tp, Delta, fval_t);
    term2 = Uxt * tail_avg * Uxt';

    M = term1 + term2;
end

function loss = loss_total(a, tp, Delta, A, Nbasis, Nt)

    sigz = [1 0; 0 -1];
    sigx = [0 1; 1 0];
    sigy = [0 -1i; 1i 0];

    tlist = linspace(0.01, tp, Nt);
    loss_matrix = zeros(4,4);
    
    parfor k = 1:Nt
        t = tlist(k);
        M = compute_M(t, a, tp, Delta, A);
        H_0 = -sigz/2;
        fval = f_s(t, a, tp);
        U2 = expm(-1i * fval/2 * [0 1; 1 0]);
        Uint = exp(-1i*H_0*t)*U2;
        Mint = Uint*M*Uint';
        Aint = Uint*A*Uint';  % A(t)
        D = kron(conj(Mint), Aint) + kron(conj(Aint), Mint) - kron(eye(2) , (Aint * Mint)) - kron((conj(Aint) * conj(Mint)) ,eye(2));

        U = (1/sqrt(2)) * [1 0 0 1; 0 1 1i 0; 0 1 -1i 0; 1 0 0 -1];
        Dp = U' * D * U;
        loss_matrix = loss_matrix + D;
    end

    loss = real(trace(loss_matrix'*loss_matrix)) / Nt;
    %loss = -real(trace(loss_matrix)) / Nt;
end
