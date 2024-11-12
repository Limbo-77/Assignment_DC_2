clear;
M_16 = 16; % 16-FSK�Ľ���
M_8 = 8;   % 8-FSK�Ľ���
M_4 = 4;   % 4-FSK�Ľ���
M_2 = 2;   % 2-FSK�Ľ���
Eb_No_dB_range = 0:15; % Eb/N0 ��Χ��0��15 dB
N = 100000; % �������

% ��ʼ���洢����
Ps_sim_values_16 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_16 = zeros(1, length(Eb_No_dB_range));
Ps_sim_values_8 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_8 = zeros(1, length(Eb_No_dB_range));
Ps_sim_values_4 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_4 = zeros(1, length(Eb_No_dB_range));
Ps_sim_values_2 = zeros(1, length(Eb_No_dB_range));
Ps_exact_values_2 = zeros(1, length(Eb_No_dB_range));

% �����������������
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

% ��ͼ
semilogy(Eb_No_dB_range, Ps_sim_values_16, 'o-', 'DisplayName', '16-FSK Simulation');
hold on;
semilogy(Eb_No_dB_range, Ps_exact_values_16, 'x-', 'DisplayName', '16-FSK Nearest Neighbor Approximation');
semilogy(Eb_No_dB_range, Ps_sim_values_8, 's-', 'DisplayName', '8-FSK Simulation');
semilogy(Eb_No_dB_range, Ps_exact_values_8, 'd-', 'DisplayName', '8-FSK Nearest Neighbor Approximation');
semilogy(Eb_No_dB_range, Ps_sim_values_4, '^-', 'DisplayName', '4-FSK Simulation');
semilogy(Eb_No_dB_range, Ps_exact_values_4, 'v-', 'DisplayName', '4-FSK Nearest Neighbor Approximation');
semilogy(Eb_No_dB_range, Ps_sim_values_2, 'p-', 'DisplayName', '2-FSK Simulation');
semilogy(Eb_No_dB_range, Ps_exact_values_2, 'h-', 'DisplayName', '2-FSK Nearest Neighbor Approximation');

% ����ͼ������
xlabel('Eb/N0 (dB)');
ylabel('Probability of Symbol Error P');
xlim([0 15]); % ���᷶Χ��0��15
ylim([1e-6 1]); % ���᷶Χ��10^(-6)��10^0
legend;
grid on;
title('16-FSK, 8-FSK, 4-FSK, and 2-FSK Symbol Error Rate vs Eb/N0');
hold off;

% ���� M-FSK ������������
function [Ps_sim, Ps_exact] = compute_MFSK_error(Eb_No_dB, M, N)
    % ����������
    constellation_point = eye(M);
    
    % �����źŹ��ʺ�������׼��
    Es = mean(abs(constellation_point(1)).^2);
    Es = Es(1);
    Eb_No_linear = 10^(Eb_No_dB/10);
    sigma = sqrt(Es/(2*Eb_No_linear*log2(M)));

    % ģ����Ŵ������
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

    % ���Ŵ�����ʵĽ��Ƽ��� (Nearest Neighbor Approximation)
    Es_No = Eb_No_linear * log2(M);
    mean_y = sqrt(2 * Es_No);
    integrand = @(y) (1 - qfunc(y)).^(M - 1) .* exp(-0.5 * (y - mean_y).^2) / sqrt(2 * pi);
    integral_value = integral(integrand, -Inf, Inf, 'ArrayValued', true, 'RelTol',1e-8, 'AbsTol',1e-12);
    Ps_exact = 1 - integral_value;
end
