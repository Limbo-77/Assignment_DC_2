clear;
M_values = [16, 8, 4, 2]; % FSK modulation orders
Eb_No_dB_range = 0:15; % Eb/N0 range from 0 to 15 dB
N = 1000000; % Simulation count

% Initialize storage for error probabilities
Ps_exact_values = cell(1, length(M_values));
Ps_NNA_values = cell(1, length(M_values));

% Loop over each M value
for m = 1:length(M_values)
    M = M_values(m);
    Ps_exact_values{m} = zeros(1, length(Eb_No_dB_range));
    Ps_NNA_values{m} = zeros(1, length(Eb_No_dB_range));
    
    for idx = 1:length(Eb_No_dB_range)
        Eb_No_dB = Eb_No_dB_range(idx);
        
        % Calculate exact performance using numerical integration
        Ps_exact_values{m}(idx) = exact_MFSK_performance(Eb_No_dB, M);
        
        % Calculate Nearest Neighbor Approximation
        [~, Ps_NNA_values{m}(idx)] = compute_MFSK_error(Eb_No_dB, M, N);
    end
end

% Plotting the results
figure;
hold on;
for m = 1:length(M_values)
    M = M_values(m);
    semilogy(Eb_No_dB_range, Ps_exact_values{m}, 'DisplayName', sprintf('%d-FSK Exact Performance (Numerical Integration)', M));
    semilogy(Eb_No_dB_range, Ps_NNA_values{m}, '--', 'DisplayName', sprintf('%d-FSK Nearest Neighbor Approximation', M));
end

% Set plot properties
xlabel('Eb/N0 (dB)');
ylabel('Probability of Symbol Error P');
xlim([0 15]);
ylim([1e-6 1]);
legend show;
grid on;
title('Comparison of Exact Performance vs Nearest Neighbor Approximation for 16-FSK, 8-FSK, 4-FSK, and 2-FSK');
hold off;

% Exact MFSK performance using numerical integration
function Ps_exact = exact_MFSK_performance(Eb_No_dB, M)
    Eb_No = 10^(Eb_No_dB / 10);
    Es_No = Eb_No * log2(M); % Calculate Es/No
    mean_y = sqrt(2 * Es_No);
    
    % Define integrand
    integrand = @(y) (1 - qfunc(y)).^(M - 1) .* exp(-0.5 * (y - mean_y).^2) / sqrt(2 * pi);
    
    % Perform numerical integration
    Ps_exact = 1 - integral(integrand, -Inf, Inf, 'ArrayValued', true, 'RelTol',1e-8, 'AbsTol',1e-12);
end
