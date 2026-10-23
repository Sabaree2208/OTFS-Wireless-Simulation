clear; clc; close all;
rng(42); % Locks random seed for consistent results

%% ============ USER SETTINGS ============
N_grid         = 16;          % Increased Doppler dimension for better resolution
M_grid         = 16;          % Increased Delay dimension to avoid path overlapping
modOrder       = 4;           % QAM modulation order (4-QAM for a clean waterfall)
cp_len         = N_grid/4;    % Cyclic prefix length
numSymbols     = N_grid * M_grid;
num_bits_per_frame = numSymbols * log2(modOrder);
num_frames     = 500;         % Number of frames for statistical averaging
snr_db_range   = -6:2:14;     % SNR range for simulation
nakagami_m     = 3.0;         % Nakagami shape parameter for milder fading
%% ======================================

% Nakagami envelope generator
if exist('gamrnd','file') == 2
    nak = @(m,Om,n) sqrt(gamrnd(m, Om/m, n, 1));
else
    nak = @(m,Om,n) sqrt((Om/(2*m)) * sum(randn(n, round(2*m)).^2, 2));
end

% 3GPP Extended Vehicular A (EVA) Channel Parameters
eva_delays_ns = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510];
eva_powers_dB = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
num_paths     = length(eva_delays_ns);

eva_power_linear = 10.^(eva_powers_dB / 10);
eva_power_normalized = eva_power_linear / sum(eva_power_linear); 

path_delays = mod(round(eva_delays_ns / max(eva_delays_ns) * (M_grid-1)), M_grid); 
path_dopplers = randi([0, N_grid-1], 1, num_paths); 
cconv2 = @(A, B) ifft2(fft2(A) .* fft2(B));

channel_types = {'3GPP', 'Nakagami'};
ber_results = zeros(2, length(snr_db_range));

%% ========================================================================
%  SIMULATION LOOP (3GPP vs Nakagami)
% =========================================================================
for ch_idx = 1:2
    disp(['--- Simulating OTFS over ' channel_types{ch_idx} ' Channel ---']);
    
    path_gains = zeros(1, num_paths);
    for p = 1:num_paths
        Om = eva_power_normalized(p);
        if ch_idx == 1
            m_val = 1.0; % Rayleigh equivalent for 3GPP baseline
        else
            m_val = nakagami_m;
        end
        envelope = nak(m_val, Om, 1);
        phase = 2 * pi * rand();
        path_gains(p) = envelope * exp(1i * phase);
    end
    
    h_dd = zeros(M_grid, N_grid);
    for i = 1:num_paths
        h_dd(path_delays(i)+1, path_dopplers(i)+1) = h_dd(path_delays(i)+1, path_dopplers(i)+1) + path_gains(i); 
    end
    
    ber_values = zeros(size(snr_db_range));
    
    for snr_idx = 1:length(snr_db_range)
        current_snr_db = snr_db_range(snr_idx);
        total_bit_errors = 0;
        
        for frame = 1:num_frames
            % --- Transmitter (Tx) ---
            message_signal_binary = randi([0 1], num_bits_per_frame, 1);
            qam_symbols = qammod(message_signal_binary, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);
            x_dd = reshape(qam_symbols, M_grid, N_grid);
            X_tf = fft(ifft(x_dd, M_grid, 1), N_grid, 2);
            x_t_matrix = ifft(X_tf, N_grid, 2);
            x_t_with_cp = [x_t_matrix(:, N_grid-cp_len+1:N_grid), x_t_matrix];
            x_t = x_t_with_cp.';
            x_t = x_t(:);
            
            % --- Channel ---
            y_dd_clean = cconv2(h_dd, x_dd);
            Y_tf_clean = fft(ifft(y_dd_clean, M_grid, 1), N_grid, 2);
            y_t_matrix = ifft(Y_tf_clean, N_grid, 2);
            y_t_with_cp = [y_t_matrix(:, N_grid-cp_len+1:N_grid), y_t_matrix];
            y_t_clean = y_t_with_cp.';
            y_t_clean = y_t_clean(:);
            y_t_noisy = awgn(y_t_clean, current_snr_db, 'measured');
        
            % --- Receiver (Rx) ---
            y_t_reshaped = reshape(y_t_noisy, N_grid + cp_len, M_grid).';
            y_t_no_cp = y_t_reshaped(:, cp_len+1:end);
            Y_tf_received = fft(y_t_no_cp, N_grid, 2);
            y_dd_received = fft(ifft(Y_tf_received, N_grid, 2), M_grid, 1);
            
            % --- 2D MMSE Equalizer ---
            snr_linear = 10^(current_snr_db / 10);
            noise_var = 1 / snr_linear; 
            
            H_f = fft2(h_dd);
            Y_f = fft2(y_dd_received);
            X_f_est = (conj(H_f) .* Y_f) ./ (abs(H_f).^2 + noise_var);
            
            y_dd_equalized = ifft2(X_f_est);
            
            % --- Demodulation ---
            qam_symbols_received = reshape(y_dd_equalized, numSymbols, 1);
            recovered_binary_signal = qamdemod(qam_symbols_received, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
            
            total_bit_errors = total_bit_errors + biterr(message_signal_binary, recovered_binary_signal);
        end
        
        ber_values(snr_idx) = total_bit_errors / (num_frames * num_bits_per_frame);
    end
    ber_results(ch_idx, :) = ber_values;
end
disp('Simulations Complete.');

%% ========================================================================
%  PLOT BER COMPARISON (3GPP in Black, Nakagami in Red)
% =========================================================================
figure('Name', 'OTFS BER Comparison: 3GPP vs Nakagami', 'Position', [300, 300, 800, 600]);
semilogy(snr_db_range, ber_results(1, :), '-o', 'Color', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
semilogy(snr_db_range, ber_results(2, :), '-s', 'Color', 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
grid on;
title('OTFS Bit Error Rate (BER) Waterfall Comparison: 3GPP vs. Nakagami');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend(channel_types, 'Location', 'southwest');
ylim([1e-5 1]);
xlim([min(snr_db_range) max(snr_db_range)]);