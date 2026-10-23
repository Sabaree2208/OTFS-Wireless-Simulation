clear; clc; close all;

%% ---------------- 0) Create Figure and Plot Block Diagram ----------------
figure('Name', 'OTFS Full System Visualization with BER', 'Position', [50, 50, 1500, 950]);
annotation('textbox', [0.02, 0.93, 0.96, 0.05], 'String', 'OTFS Block Diagram and Step-by-Step Signal Flow', ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'EdgeColor', 'none');

% Define positions for the diagram blocks
x_pos = 0.04; y_pos = 0.84; w = 0.075; h = 0.07; g = 0.008;
blocks = {'1. Message Bits', '2. QAM Modulation', '3. DD Grid Mapper', '4. ISFFT (DD→TF)', ...
    '5. Heisenberg (TF→Time)', 'Channel', '6. Wigner (Time→TF)', '7. SFFT (TF→DD)', ...
    '8. DD Equalizer', '9. QAM Demodulator','10. Recovered Bits'};

% Loop to draw each block and the arrow that follows it
for i = 1:numel(blocks)
    annotation('rectangle', [x_pos, y_pos, w, h], 'FaceColor', [0.9, 0.95, 1.0]);
    annotation('textbox', [x_pos, y_pos, w, h], 'String', blocks{i}, 'EdgeColor', 'none', ...
        'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    if i < numel(blocks)
        annotation('arrow', [x_pos + w, x_pos + w + g], [y_pos + h/2, y_pos + h/2]);
    end
    x_pos = x_pos + w + g;
end

%% --- 1. System and Channel Definition ---
% System Parameters
N = 8; M = 8; modOrder = 4; cp_len = N/4;
numSymbols = N * M;
num_bits_per_frame = numSymbols * log2(modOrder);
num_frames = 5000;   9   % Number of frames for statistical averaging
snr_db_range = -6:2:14; % SNR range for simulation

% --- 3GPP Extended Vehicular A (EVA) Channel Parameters ---
eva_delays_ns = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510];
eva_powers_dB = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];

% Convert dB power to linear amplitudes and apply Rayleigh fading
eva_power_linear = 10.^(eva_powers_dB / 10);
eva_power_normalized = eva_power_linear / sum(eva_power_linear); % Normalize total power
rayleigh_fading = (randn(1, 9) + 1i * randn(1, 9)) / sqrt(2);
path_gains = sqrt(eva_power_normalized) .* rayleigh_fading;

% Map the 9 paths to your specific Delay-Doppler Grid (M=8, N=8)
% (Scales the 2510ns max delay to fit the M=8 grid and randomizes Doppler)
path_delays = mod(round(eva_delays_ns / max(eva_delays_ns) * (M-1)), M); 
path_dopplers = randi([0, N-1], 1, 9); 

% Populate the Delay-Doppler channel matrix (h_dd)
h_dd = zeros(M, N);
for i = 1:length(path_gains)
    % Accumulate gains if multiple paths map to the same Delay-Doppler bin
    h_dd(path_delays(i)+1, path_dopplers(i)+1) = h_dd(path_delays(i)+1, path_dopplers(i)+1) + path_gains(i); 
end

% Helper function for 2D circular convolution
cconv2 = @(A, B) ifft2(fft2(A) .* fft2(B));

%% ========================================================================
%  PART 2: Main Simulation Loop for BER Calculation
% =========================================================================
disp('--- Simulating Transmitter and Receiver for BER ---');
ber_values = zeros(size(snr_db_range));

for snr_idx = 1:length(snr_db_range)
    current_snr_db = snr_db_range(snr_idx);
    total_bit_errors = 0;
    
    for frame = 1:num_frames
        % --- Transmitter (Tx) ---
        message_signal_binary = randi([0 1], num_bits_per_frame, 1);
        qam_symbols = qammod(message_signal_binary, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);
        x_dd = reshape(qam_symbols, M, N);
        X_tf = fft(ifft(x_dd, M, 1), N, 2);
        x_t_matrix = ifft(X_tf, N, 2);
        x_t_with_cp = [x_t_matrix(:, N-cp_len+1:N), x_t_matrix];
        x_t = x_t_with_cp.';
        x_t = x_t(:);
        
        % --- Channel ---
        y_dd_clean = cconv2(h_dd, x_dd);
        Y_tf_clean = fft(ifft(y_dd_clean, M, 1), N, 2);
        y_t_matrix = ifft(Y_tf_clean, N, 2);
        y_t_with_cp = [y_t_matrix(:, N-cp_len+1:N), y_t_matrix];
        y_t_clean = y_t_with_cp.';
        y_t_clean = y_t_clean(:);
        y_t_noisy = awgn(y_t_clean, current_snr_db, 'measured');
    
       
        % --- Receiver (Rx) ---
        y_t_reshaped = reshape(y_t_noisy, N + cp_len, M).';
        y_t_no_cp = y_t_reshaped(:, cp_len+1:end);
        Y_tf_received = fft(y_t_no_cp, N, 2);
        y_dd_received = fft(ifft(Y_tf_received, N, 2), M, 1);
        
        % --- 2D MMSE Equalizer ---
        % Convert SNR to linear scale to estimate noise variance
        snr_linear = 10^(current_snr_db / 10);
        noise_var = 1 / snr_linear; 
        
        % Perform 2D deconvolution in the frequency domain
        H_f = fft2(h_dd);
        Y_f = fft2(y_dd_received);
        X_f_est = (conj(H_f) .* Y_f) ./ (abs(H_f).^2 + noise_var);
        
        % Return to Delay-Doppler domain
        y_dd_equalized = ifft2(X_f_est);
        
        % --- Demodulation ---
        qam_symbols_received = reshape(y_dd_equalized, numSymbols, 1);
        recovered_binary_signal = qamdemod(qam_symbols_received, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        qam_symbols_received = reshape(y_dd_equalized, numSymbols, 1);
        recovered_binary_signal = qamdemod(qam_symbols_received, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);
        
        % Accumulate errors for this frame
        total_bit_errors = total_bit_errors + biterr(message_signal_binary, recovered_binary_signal);
    end
    
    % Calculate average BER for this SNR point
    ber_values(snr_idx) = total_bit_errors / (num_frames * num_bits_per_frame);
    fprintf('SNR: %d dB, Total Bit Errors: %d, BER: %e\n', current_snr_db, total_bit_errors, ber_values(snr_idx));
end
disp('BER Simulation Complete.');

%% ========================================================================
%  PART 3: Visualization of the FINAL FRAME from the FINAL SNR run
% =========================================================================
% The variables from the last iteration (x_dd, X_tf, etc.) are used for plotting.

% --- Tx Plots (1-5) ---
custom_subplot(4, 3, 1); stem(message_signal_binary, 'filled'); grid on; title('1. Original Message Signal (Last Frame)'); xlabel('Bit Index'); ylabel('Value'); xlim([0 num_bits_per_frame]); ylim([-0.2 1.2]);
custom_subplot(4, 3, 2); plot(real(qam_symbols), imag(qam_symbols), 'o', 'MarkerFaceColor', 'b'); grid on; title('2. QAM Modulated Symbols'); xlabel('In-Phase'); ylabel('Quadrature'); axis([-2 2 -2 2]); axis square;
custom_subplot(4, 3, 3); surf(0:N-1, 0:M-1, abs(x_dd)); title({'3. Signal on Delay-Doppler Grid', '(x[k,l])'}); xlabel('Doppler Index (k)'); ylabel('Delay Index (l)'); zlabel('Magnitude'); colormap('parula'); axis square; shading interp;
custom_subplot(4, 3, 4); surf(0:N-1, 0:M-1, abs(X_tf)); title({'4. Signal in Time-Frequency Domain', '(X[n,m])'}); xlabel('Frequency Index (n)'); ylabel('Time Index (m)'); zlabel('Magnitude'); axis square; shading interp;
custom_subplot(4, 3, 5); plot(abs(x_t)); grid on; title({'5. Transmitted Time-Domain Signal', '(x(t))'}); xlabel('Time Samples'); ylabel('Magnitude'); xlim([0 length(x_t)]);

% --- Rx Plots (6-11) ---
custom_subplot(4, 3, 6); plot(abs(y_t_noisy)); grid on; title({'6. Received Time-Domain Signal (Last SNR)', '(y(t))'}); xlabel('Time Samples'); ylabel('Magnitude'); xlim([0 length(y_t_noisy)]);
custom_subplot(4, 3, 7); surf(0:N-1, 0:M-1, abs(Y_tf_received)); title({'7. Received Signal in TF Domain', '(Y[n,m])'}); xlabel('Frequency Index (n)'); ylabel('Time Index (m)'); zlabel('Magnitude'); axis square; shading interp;
custom_subplot(4, 3, 8); surf(0:N-1, 0:M-1, abs(y_dd_received)); title({'8. Received Signal in DD Domain', '(y[k,l])'}); xlabel('Doppler Index (k)'); ylabel('Delay Index (l)'); zlabel('Magnitude'); axis square; shading interp;
custom_subplot(4, 3, 9); plot(real(qam_symbols_received), imag(qam_symbols_received), 'x', 'Color', [0.8500 0.3250 0.0980]); hold on; plot(real(qam_symbols), imag(qam_symbols), 'o', 'MarkerEdgeColor', 'b'); grid on; title('9. Equalized QAM Symbols'); xlabel('In-Phase'); ylabel('Quadrature'); axis([-2 2 -2 2]); axis square; legend('Received', 'Original', 'Location', 'southeast');
custom_subplot(4, 3, 10); stem(recovered_binary_signal, 'filled', 'Color', [0.4660 0.6740 0.1880]); grid on; title('10. Recovered Message Signal (Binary)'); xlabel('Bit Index'); ylabel('Value'); xlim([0 num_bits_per_frame]); ylim([-0.2 1.2]);
custom_subplot(4, 3, 11); stem(message_signal_binary, 'b', 'DisplayName', 'Original'); hold on; stem(recovered_binary_signal, 'r--', 'DisplayName', 'Recovered'); grid on; title('11. Final Comparison (Last Frame)'); xlabel('Bit Index'); ylabel('Value'); xlim([0 num_bits_per_frame]); ylim([-0.2 1.2]); legend('show');

% --- Final BER Plot (12) ---
custom_subplot(4, 3, 12);
semilogy(snr_db_range, ber_values, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; title('12. Bit Error Rate (BER) vs. SNR');
xlabel('SNR (dB)'); ylabel('Bit Error Rate');
ylim([1e-5 1]); xlim([min(snr_db_range) max(snr_db_range)]);

%% ========================================================================
%  HELPER FUNCTION FOR CUSTOM SUBPLOTTING
% =========================================================================
function custom_subplot(rows, cols, index)
    % This function manually calculates the position for each subplot.
    top_margin=0.22; bottom_margin=0.05; left_margin=0.05; right_margin=0.05;
    horz_padding=0.06; vert_padding=0.08;
    total_plot_h = 1.0 - top_margin - bottom_margin;
    total_plot_w = 1.0 - left_margin - right_margin;
    plot_w = (total_plot_w - (cols-1) * horz_padding) / cols;
    plot_h = (total_plot_h - (rows-1) * vert_padding) / rows;
    current_col = mod(index - 1, cols);
    current_row = floor((index - 1) / cols);
    left = left_margin + current_col * (plot_w + horz_padding);
    bottom = 1.0 - top_margin - plot_h - current_row * (plot_h + vert_padding);
    axes('Position', [left, bottom, plot_w, plot_h]);
end