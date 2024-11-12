clear;
M_16 = 16; % 16-FSK的阶数
M_8 = 8;   % 8-FSK的阶数
M_4 = 4;   % 4-FSK的阶数
M_2 = 2;   % 2-FSK的阶数
Eb_No_dB_range = 0:15; % Eb/N0 范围从0到15 dB
N = 100000; % 仿真次数

% 初始化存储变量
Ps_sim_values_16 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_16 = zeros(1, length(Eb_No_dB_range));
Ps_sim_values_8 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_8 = zeros(1, length(Eb_No_dB_range));
Ps_sim_values_4 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_4 = zeros(1, length(Eb_No_dB_range));
Ps_sim_values_2 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_2 = zeros(1, length(Eb_No_dB_range));

% 计算各阶数的误差概率
for idx = 1:length(Eb_No_dB_range)
    Eb_No_dB = Eb_No_dB_range(idx);
    
    % 16-FSK
    [Ps_sim_values_16(idx), Ps_exact_values_16(idx)] = compute_MFSK_error(Eb_No_dB, M_16, N);
    
    % 8-FSK
    [Ps_sim_values_8(idx), Ps_exact_values_8(idx)] = compute_MFSK_error(Eb_No_dB, M_8, N);
    
    % 4-FSK
    [Ps_sim_values_4(idx), Ps_exact_values_4(idx)] = compute_MFSK_error(Eb_No_dB, M_4, N);
    
    % 2-FSK
    [Ps_sim_values_2(idx), Ps_exact_values_2(idx)] = compute_MFSK_error(Eb_No_dB, M_2, N);
end

% 绘图
semilogy(Eb_No_dB_range, Ps_sim_values_16, 'o-', 'DisplayName', '16-FSK Simulation');
hold on;
semilogy(Eb_No_dB_range, Ps_exact_values_16, 'x-', 'DisplayName', '16-FSK Nearest Neighbor Approximation');
semilogy(Eb_No_dB_range, Ps_sim_values_8, 's-', 'DisplayName', '8-FSK Simulation');
semilogy(Eb_No_dB_range, Ps_exact_values_8, 'd-', 'DisplayName', '8-FSK Nearest Neighbor Approximation');
semilogy(Eb_No_dB_range, Ps_sim_values_4, '^-', 'DisplayName', '4-FSK Simulation');
semilogy(Eb_No_dB_range, Ps_exact_values_4, 'v-', 'DisplayName', '4-FSK Nearest Neighbor Approximation');
semilogy(Eb_No_dB_range, Ps_sim_values_2, 'p-', 'DisplayName', '2-FSK Simulation');
semilogy(Eb_No_dB_range, Ps_exact_values_2, 'h-', 'DisplayName', '2-FSK Nearest Neighbor Approximation');

% 设置图表属性
xlabel('Eb/N0 (dB)');
ylabel('Probability of Symbol Error P');
xlim([0 15]); % 横轴范围从0到15
ylim([1e-6 1]); % 纵轴范围从10^(-6)到10^0
legend;
grid on;
title('16-FSK, 8-FSK, 4-FSK, and 2-FSK Symbol Error Rate vs Eb/N0');
hold off;

% 计算 M-FSK 误差的内联函数
function [Ps_sim, Ps_exact] = compute_MFSK_error(Eb_No_dB, M, N)
    % 生成星座点
    constellation_point = eye(M);
    
    % 计算信号功率和噪声标准差
    Es = mean(abs(constellation_point(1)).^2);
    Es = Es(1);
    Eb_No_linear = 10^(Eb_No_dB/10);
    sigma = sqrt(Es/(2*Eb_No_linear*log2(M)));

    % 模拟符号错误概率
    number_of_symbol_error = 0;
    for i = 1:N
        transmitted_index = 1 + floor(rand * M);
        BB_signal = constellation_point(transmitted_index,:);
        received_signal = BB_signal + (randn(1, M) .* sigma);
        [~, decision_index] = max(received_signal);
        if decision_index ~= transmitted_index
            number_of_symbol_error = number_of_symbol_error + 1;
        end
    end
    Ps_sim = number_of_symbol_error / N;

    % 符号错误概率的近似计算 (Nearest Neighbor Approximation)
    Es_No = Eb_No_linear * log2(M);
    mean_y = sqrt(2 * Es_No);
    integrand = @(y) (1 - qfunc(y)).^(M - 1) .* exp(-0.5 * (y - mean_y).^2) / sqrt(2 * pi);
    integral_value = integral(integrand, -Inf, Inf, 'ArrayValued', true, 'RelTol',1e-8, 'AbsTol',1e-12);
    Ps_exact = 1 - integral_value;
end
