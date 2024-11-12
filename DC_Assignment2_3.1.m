clear;
M16 = 16;
M4 = 4;
d = 1;

% Define 16-QAM constellation
constellation_point_16QAM = [-3*d + 1j*3*d  -d + 1j*3*d   d + 1j*3*d   3*d + 1j*3*d ...
                             -3*d + 1j*d    -d + 1j*d     d + 1j*d     3*d + 1j*d ...
                             -3*d - 1j*d    -d - 1j*d     d - 1j*d     3*d - 1j*d ...
                             -3*d - 1j*3*d  -d - 1j*3*d   d - 1j*3*d   3*d - 1j*3*d];
Es16 = mean(abs(constellation_point_16QAM).^2);

% Define 4-QAM constellation
constellation_point_4QAM = [-d + 1j*d, d + 1j*d, -d - 1j*d, d - 1j*d];
Es4 = mean(abs(constellation_point_4QAM).^2);

% Define 16-PSK constellation
constellation_point_16PSK = exp(1j * (0:M16-1) * 2 * pi / M16);
Es16PSK = mean(abs(constellation_point_16PSK).^2);

Eb_No_db = 0:20; % Eb/N0 范围从0到20 dB
N = 10000; % 增加符号数以提高准确性

% Initialize arrays for storing symbol error rates
Ps_sim_16QAM = zeros(1, length(Eb_No_db));
Ps_exact_16QAM = zeros(1, length(Eb_No_db));
Ps_sim_4QAM = zeros(1, length(Eb_No_db));
Ps_exact_4QAM = zeros(1, length(Eb_No_db));
Ps_sim_16PSK = zeros(1, length(Eb_No_db));
Ps_exact_16PSK = zeros(1, length(Eb_No_db));

% Simulation and Nearest Neighbor Approximation for 16-QAM
for idx = 1:length(Eb_No_db)
    number_of_symbol_error_16QAM = 0;
    Eb_No_linear = 10^(Eb_No_db(idx)/10);
    sgma16 = sqrt(Es16 / (2 * Eb_No_linear * log2(M16)));
    
    % 16-QAM Simulation
    for i = 1:N
        transmitted_index = 1 + floor(rand * M16);
        BB_signal = constellation_point_16QAM(transmitted_index);
        
        received_signal = BB_signal + (randn * sgma16 + 1j * randn * sgma16);
        
        % Minimum Distance Detection for 16-QAM
        distances = abs(received_signal - constellation_point_16QAM).^2;
        [~, decision_index] = min(distances);
        
        if decision_index ~= transmitted_index
            number_of_symbol_error_16QAM = number_of_symbol_error_16QAM + 1;
        end
    end
    Ps_sim_16QAM(idx) = number_of_symbol_error_16QAM / N;
    
    % Nearest Neighbor Approximation for 16-QAM
    Ps_exact_16QAM(idx) = 1 - (1 - 2 * (sqrt(M16) - 1) / sqrt(M16) * qfunc(sqrt(3 * Eb_No_linear * log2(M16) / (M16 - 1))))^2;
end

% Simulation and Nearest Neighbor Approximation for 4-QAM
for idx = 1:length(Eb_No_db)
    number_of_symbol_error_4QAM = 0;
    Eb_No_linear = 10^(Eb_No_db(idx)/10);
    sgma4 = sqrt(Es4 / (2 * Eb_No_linear * log2(M4)));
    
    % 4-QAM Simulation
    for i = 1:N
        transmitted_index = 1 + floor(rand * M4);
        BB_signal = constellation_point_4QAM(transmitted_index);
        
        received_signal = BB_signal + (randn * sgma4 + 1j * randn * sgma4);
        
        % Minimum Distance Detection for 4-QAM
        distances = abs(received_signal - constellation_point_4QAM).^2;
        [~, decision_index] = min(distances);
        
        if decision_index ~= transmitted_index
            number_of_symbol_error_4QAM = number_of_symbol_error_4QAM + 1;
        end
    end
    Ps_sim_4QAM(idx) = number_of_symbol_error_4QAM / N;
    
    % Nearest Neighbor Approximation for 4-QAM
    Ps_exact_4QAM(idx) = 2 * qfunc(sqrt(2 * Eb_No_linear));
end

% Simulation and Nearest Neighbor Approximation for 16-PSK
for idx = 1:length(Eb_No_db)
    number_of_symbol_error_16PSK = 0;
    Eb_No_linear = 10^(Eb_No_db(idx)/10);
    sgma16PSK = sqrt(Es16PSK / (2 * Eb_No_linear * log2(M16)));
    
    % 16-PSK Simulation
    for i = 1:N
        transmitted_index = 1 + floor(rand * M16);
        BB_signal = constellation_point_16PSK(transmitted_index);
        
        received_signal = BB_signal + (randn * sgma16PSK + 1j * randn * sgma16PSK);
        
        % Minimum Distance Detection for 16-PSK
        distances = abs(received_signal - constellation_point_16PSK).^2;
        [~, decision_index] = min(distances);
        
        if decision_index ~= transmitted_index
            number_of_symbol_error_16PSK = number_of_symbol_error_16PSK + 1;
        end
    end
    Ps_sim_16PSK(idx) = number_of_symbol_error_16PSK / N;
    
    % Nearest Neighbor Approximation for 16-PSK
    Ps_exact_16PSK(idx) = 2 * qfunc(sqrt(2 * Eb_No_linear * log2(M16)) * sin(pi / M16));

end

% Plotting
semilogy(Eb_No_db, Ps_sim_16QAM, 'o-', 'DisplayName', '16-QAM Simulation');
hold on;
semilogy(Eb_No_db, Ps_exact_16QAM, 'x-', 'DisplayName', '16-QAM Nearest Neighbor Approximation');
semilogy(Eb_No_db, Ps_sim_4QAM, 's-', 'DisplayName', '4-QAM Simulation');
semilogy(Eb_No_db, Ps_exact_4QAM, 'd-', 'DisplayName', '4-QAM Nearest Neighbor Approximation');
semilogy(Eb_No_db, Ps_sim_16PSK, '^-', 'DisplayName', '16-PSK Simulation');
semilogy(Eb_No_db, Ps_exact_16PSK, 'p-', 'DisplayName', '16-PSK Nearest Neighbor Approximation');
xlabel('Eb/N0 (dB)');
ylabel('Probability of Symbol Error P');
xlim([0 20]); % 设置横轴范围从0到20
ylim([1e-8 1]); % 设置纵轴范围从10^(-8)到10^0
legend;
grid on;
title('Symbol Error Rate vs Eb/N0 for 16-QAM, 4-QAM, and 16-PSK');
hold off;

