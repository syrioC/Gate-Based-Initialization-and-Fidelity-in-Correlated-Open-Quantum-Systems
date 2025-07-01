
N = 48;
fid_data = zeros(N, 3);  % columns: xi, fid_opt, fid_base
xi_vals = linspace(0, 10, N);
delete(gcp('nocreate'));
parpool('local', 8);
parfor k = 1:N
    xi = xi_vals(k);
    k
    [rho_opt, rho_base] = simulate_fidelity(1e-05, xi);
    [fid_opt,~,~] = bloch_coordinates_from_rho(rho_opt);
    [fid_base,~,~] = bloch_coordinates_from_rho(rho_base);
    fid_data(k, :) = [xi, fid_opt, fid_base];
end

writematrix(fid_data, 'fidelity_vs_xi_version3.txt', 'Delimiter', '\t');
figure(1);
plot(xi_vals, fid_data(:,2),'b-'); hold on
plot(xi_vals, fid_data(:,3),'r-'); hold on
delete(gcp('nocreate'));

function [fid, theta, phi] = bloch_coordinates_from_rho(rho)
    % Pauli matrices
    sigx = [0 1; 1 0];
    sigy = [0 -1i; 1i 0];
    sigz = [1 0; 0 -1];

    % Bloch vector components
    x = real(trace(rho * sigx));
    y = real(trace(rho * sigy));
    z = real(trace(rho * sigz));

    % Radius (purity measure)
    r = sqrt(x^2 + y^2 + z^2);

    % Avoid division by zero
    if r > 1e-10
        theta = acos(z / r);
        phi = atan2(y, x);
    else
        theta = 0;
        phi = 0;
    end
    x = sin(theta) * cos(phi);
    y = sin(theta) * sin(phi);
    z = cos(theta);
            
    % Construct the pure state density matrix for this Bloch vector
    bloch_vector = [x; y; z];
    rho_ideal = 0.5 * (eye(2) + bloch_vector(1) * [0 1; 1 0] + ...
                                bloch_vector(2) * [0 -1i; 1i 0] + ...
                                bloch_vector(3) * [1 0; 0 -1]);
    fid = 1 - 0.5 * sum(svd(rho - rho_ideal));
end
