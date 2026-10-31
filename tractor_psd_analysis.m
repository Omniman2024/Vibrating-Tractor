function tractor_psd_analysis()
    % SYSTEM PARAMETERS
    mc = 2700;         
    Jc = 3600;          
    kf = 714400;        
    kr = 718600;       
    cf = 3600;         
    cr = 5100;          
    lf = 1.35;          
    lr = 0.85;          
    v = 6 / 3.6;        
    
    % SOLVER SETUP
    x0 = [0; 0; 0; 0];  
    tspan = 0:0.01:10; 

    % NUMERICAL INTEGRATION (Runge-Kutta 4th/5th order)
    [t, X] = ode45(@(t,x) state_space_matrix(t, x, mc, Jc, kf, kr, cf, cr, lf, lr, v), tspan, x0);

    % 2. ACCELERATION RECOVERY
    zc_ddot = zeros(length(t), 1);
    phi_ddot = zeros(length(t), 1);
    
    for i = 1:length(t)
        dx = state_space_matrix(t(i), X(i,:)', mc, Jc, kf, kr, cf, cr, lf, lr, v);
        zc_ddot(i) = dx(2);    
        phi_ddot(i) = dx(4);  
    end

    % 3. CORE MATLAB FFT TO CALCULATE PSD (No Toolbox Required)
    fs = 100; % Sampling frequency (1 / 0.01 dt)
    N = length(t); % Number of samples
    
    % Vertical PSD Math
    Z_fft = fft(zc_ddot);
    Z_fft = Z_fft(1:floor(N/2)+1); % Take positive frequencies
    P_z = (1/(fs*N)) * abs(Z_fft).^2; % Calculate power
    P_z(2:end-1) = 2 * P_z(2:end-1); % Conserve energy for single-sided spectrum
    
    % Pitch PSD Math
    Phi_fft = fft(phi_ddot);
    Phi_fft = Phi_fft(1:floor(N/2)+1);
    P_phi = (1/(fs*N)) * abs(Phi_fft).^2;
    P_phi(2:end-1) = 2 * P_phi(2:end-1);
    
    % Frequency Vector
    f_psd = fs * (0:floor(N/2)) / N;

    % 4. PLOTTING THE PSD GRAPHS (With Forced Black Text)
    figure('Name', 'Power Spectral Density (PSD) Analysis', 'Color', 'w', 'Position', [150, 150, 800, 600]);
    
    ax1 = subplot(2,1,1);
    plot(f_psd, P_z, 'k', 'LineWidth', 1.5);
    xlim([0 15]); 
    xline(ax1, 3.2, 'k:', 'LineWidth', 2); 
    xline(ax1, 4.0, 'k:', 'LineWidth', 2);
    title('Vertical Vibration Power Spectral Density', 'FontWeight', 'bold', 'Color', 'k', 'FontSize', 20);
    xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'Color', 'k');
    ylabel('PSD ((m/s^2)^2/Hz)', 'FontWeight', 'bold', 'Color', 'k'); 
    % Force axes text and ticks to be black
    set(ax1, 'Color', [0.95 0.95 0.95], 'XColor', 'k', 'YColor', 'k'); 
    grid on;
    
    ax2 = subplot(2,1,2);
    plot(f_psd, P_phi, 'k', 'LineWidth', 1.5);
    xlim([0 15]);
    xline(ax2, 3.2, 'k:', 'LineWidth', 2); 
    xline(ax2, 4.0, 'k:', 'LineWidth', 2);
    title('Pitch Vibration Power Spectral Density', 'FontWeight', 'bold', 'Color', 'k', 'FontSize', 20);
    ylabel('PSD ((rad/s^2)^2/Hz)', 'FontWeight', 'bold', 'Color', 'k'); 
    xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'Color', 'k'); 
    % Force axes text and ticks to be black
    set(ax2, 'Color', [0.95 0.95 0.95], 'XColor', 'k', 'YColor', 'k'); 
    grid on;
end

% MATHEMATICAL ENGINE
function dx = state_space_matrix(t, x, mc, Jc, kf, kr, cf, cr, lf, lr, v)
    x1 = x(1); 
    x2 = x(2); 
    x3 = x(3); 
    x4 = x(4); 
    delay = (lf + lr) / v; 
    f_start = 0;
    f_end = 10;
    t_total = 10;
    freq = f_start + ((f_end - f_start) / t_total) * t;
    amp = 0.02; 
    hf = amp * sin(2*pi*freq*t); 
    freq_rate = (f_end - f_start) / t_total;
    hf_dot = amp * 2*pi * (freq + freq_rate * t) * cos(2*pi*freq*t);
    if t >= delay
        freq_del = f_start + ((f_end - f_start) / t_total) * (t-delay);
        hr = amp * sin(2*pi*freq_del*(t-delay));
        hr_dot = amp * 2*pi * (freq_del + freq_rate * (t-delay)) * cos(2*pi*freq_del*(t-delay));
    else
        hr = 0; 
        hr_dot = 0;
    end
    Fz = kf*hf + cf*hf_dot + kr*hr + cr*hr_dot; 
    Mphi = cr*lr*hr_dot - cf*lf*hf_dot + kr*lr*hr - kf*lf*hf;
    dx1 = x2;
    dx2 = (1/mc) * (Fz - (cf+cr)*x2 - (kf+kr)*x1 - (cr*lr - cf*lf)*x4 - (kr*lr - kf*lf)*x3);
    dx3 = x4;
    dx4 = (1/Jc) * (Mphi - (cr*lr - cf*lf)*x2 - (kr*lr - kf*lf)*x1 - (cf*lf^2 + cr*lr^2)*x4 - (kf*lf^2 + kr*lr^2)*x3);
    dx = [dx1; dx2; dx3; dx4];
end