function tractor_kinematics()
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
    x0 = [0; 0; 0; 0];  % Initial conditions
    tspan = 0:0.01:10;  % Simulate from t=0 to 10 seconds with 0.01s steps

    % NUMERICAL INTEGRATION (Runge-Kutta 4th/5th order)
    [t, X] = ode45(@(t,x) state_space_matrix(t, x, mc, Jc, kf, kr, cf, cr, lf, lr, v), tspan, x0);
    zc = X(:,1);        % Vertical displacement
    phi = X(:,3);      % Pitch angle     
    
    % PLOTTING
    fig = figure('Name', 'Kinematic Animation Export', 'Color', 'w', 'Position', [100, 100, 800, 450]);
    ax = axes('Parent', fig, 'Color', [0.9 0.95 1]); 
    set(ax, 'XColor', 'k', 'YColor', 'k');
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-2 2]); ylim(ax, [-0.5 0.5]);
    title(ax, 'Tractor Kinematic Animation', 'FontWeight', 'bold', 'Color', 'k','FontSize', 20);
    xlabel(ax, 'Longitudinal Distance (m)', 'FontWeight', 'bold', 'Color', 'k'); 
    ylabel(ax, 'Vertical Distance (m)', 'FontWeight', 'bold', 'Color', 'k');
    x_ground = linspace(-2, 2, 500);
    plot_ground = plot(ax, x_ground, -0.2 * ones(1, 500), 'k', 'LineWidth', 2);
    plot_chassis = plot(ax, 0, 0, 'b', 'LineWidth', 10);
    plot_com = plot(ax, 0, 0, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
    plot_wheels = plot(ax, 0, 0, 'ro', 'MarkerSize', 24, 'LineWidth', 5);
    txt_tracker = text(ax, -1.8, 0.35, 'Freq: 0.00 Hz', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
    
    % GIF EXPORT LOOP
    gif_filename = 'Tractor_Kinematics.gif';
    frame_skip = 5;
    for k = 1:frame_skip:length(t)
        t_curr = t(k);
        t_ground = t_curr + (x_ground - lf) / v; 
        y_ground = zeros(size(x_ground)); 
        valid = t_ground >= 0;           
        amp = 0.02;
        freq_ground = t_ground(valid); 
        y_ground(valid) = amp * sin(2 * pi * freq_ground .* t_ground(valid));
        set(plot_ground, 'YData', y_ground - 0.2);
        front_x = lf * cos(phi(k)); front_y = zc(k) + lf * sin(phi(k));
        rear_x = -lr * cos(phi(k)); rear_y = zc(k) - lr * sin(phi(k));
        set(plot_chassis, 'XData', [rear_x, front_x], 'YData', [rear_y, front_y]);
        set(plot_com, 'XData', 0, 'YData', zc(k));
        set(plot_wheels, 'XData', [rear_x, front_x], 'YData', [rear_y-0.2, front_y-0.2]);
        set(txt_tracker, 'String', sprintf('Freq: %.2f Hz', t(k)));
        drawnow limitrate;
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if k == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
        end
    end
    fprintf('Saved %s to your current folder.\n', gif_filename);
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