function tractor_simulation()
    % SYSTEM PARAMETERS
    mc = 2700;         
    Jc = 3600;          
    kf = 714400;        
    kr = 718600;       
    cf = 3600;         
    cr = 5100;          
    lf = 1.35;          
    lr = 0.85;          
    v = 6 / 3.6;        % Tractor speed converted from km/h to m/s

    % SOLVER SETUP
    x0 = [0; 0; 0; 0];  % Initial conditions
    tspan = 0:0.01:10;  % Simulate from t=0 to 10 seconds with 0.01s steps

    % NUMERICAL INTEGRATION (Runge-Kutta 4th/5th order)
    [t, X] = ode45(@(t,x) state_space_matrix(t, x, mc, Jc, kf, kr, cf, cr, lf, lr, v), tspan, x0);
    zc = X(:,1);        % Vertical displacement
    phi = X(:,3);      % Pitch angle
    zc_ddot = zeros(length(t), 1);
    phi_ddot = zeros(length(t), 1);
    for i = 1:length(t)
        dx = state_space_matrix(t(i), X(i,:)', mc, Jc, kf, kr, cf, cr, lf, lr, v);
        zc_ddot(i) = dx(2);    % Vertical acceleration
        phi_ddot(i) = dx(4);  % Pitch acceleration
    end
    rms_vertical_accel = rms(zc_ddot);
    rms_pitch_accel = rms(phi_ddot);
    
    fprintf('\n--- SWEPT-SINE SIMULATION RESULTS ---\n');
    fprintf('RMS Vertical Acceleration: %.2f m/s^2\n', rms_vertical_accel);
    fprintf('RMS Pitch Acceleration: %.2f rad/s^2\n\n', rms_pitch_accel);

    % PLOTTING
    figure('Name', 'Tractor 2-DOF Swept-Sine Response', 'Color', 'w');
    ax1 = subplot(2,1,1);
    plot(t, zc, 'b', 'LineWidth', 1.2);
    xline(ax1, 3.2, 'k:', 'LineWidth', 2); 
    xline(ax1, 4.0, 'k:', 'LineWidth', 2);
    ylabel('Vertical Bounce z_c (m)', 'FontWeight', 'bold', 'Color', 'k');
    xlabel('Excitation Frequency (Hz)', 'FontWeight', 'bold', 'Color', 'k');
    title('Frequency Sweep Response (0 Hz to 10 Hz)', 'FontWeight', 'bold', 'Color', 'k', 'FontSize', 20);
    set(ax1, 'Color', [0.95 0.95 0.95], 'XColor', 'k', 'YColor', 'k');
    grid on;
    ax2 = subplot(2,1,2);
    plot(t, phi, 'r', 'LineWidth', 1.2);
    xline(ax2, 3.2, 'k:', 'LineWidth', 2); 
    xline(ax2, 4.0, 'k:', 'LineWidth', 2);
    ylabel('Pitch Angle \phi_c (rad)', 'FontWeight', 'bold', 'Color', 'k');
    xlabel('Excitation Frequency (Hz)', 'FontWeight', 'bold', 'Color', 'k');
    set(ax2, 'Color', [0.95 0.95 0.95], 'XColor', 'k', 'YColor', 'k');
    grid on;
end

% MATHEMATICAL ENGINE
function dx = state_space_matrix(t, x, mc, Jc, kf, kr, cf, cr, lf, lr, v)
    x1 = x(1); % zc
    x2 = x(2); % zc_dot
    x3 = x(3); % phi
    x4 = x(4); % phi_dot
    delay = (lf + lr) / v; % Time delay for rear wheel hitting the same bump
    
    % Swept-Sine (Chirp) Road Excitation
    f_start = 0;
    f_end = 10;
    t_total = 10;
    freq = f_start + ((f_end - f_start) / t_total) * t;
    amp = 0.02; % Height of road bump = 2cm (arbitrary)
    hf = amp * sin(2*pi*freq*t); % front-wheel excitation
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

    Fz = kf*hf + cf*hf_dot + kr*hr + cr*hr_dot; % vertical force
    Mphi = cr*lr*hr_dot - cf*lf*hf_dot + kr*lr*hr - kf*lf*hf; % pitch moment

    % The 4 first-order differential equations
    dx1 = x2;
    dx2 = (1/mc) * (Fz - (cf+cr)*x2 - (kf+kr)*x1 - (cr*lr - cf*lf)*x4 - (kr*lr - kf*lf)*x3);
    dx3 = x4;
    dx4 = (1/Jc) * (Mphi - (cr*lr - cf*lf)*x2 - (kr*lr - kf*lf)*x1 - (cf*lf^2 + cr*lr^2)*x4 - (kf*lf^2 + kr*lr^2)*x3);

    dx = [dx1; dx2; dx3; dx4];
end