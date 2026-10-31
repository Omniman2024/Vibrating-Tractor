function tractor_phase_portraits()
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
    zc = X(:,1); % vertical displacement
    zc_dot = X(:,2); % vertical velocity
    phi = X(:,3); % pitch angle
    phi_dot = X(:,4); % pitch angular velocity     
    
    % PLOTTING
    fig = figure('Name', 'Phase Portrait Export', 'Color', 'w', 'Position', [200, 200, 1000, 450]);
    ax1 = subplot(1,2,1);
    set(ax1, 'Color', [0.1 0.1 0.1], 'XColor', 'k', 'YColor', 'k');
    title(ax1, 'Vertical Phase Portrait', 'FontWeight', 'bold', 'Color', 'k');
    xlabel(ax1, 'Displacement (m)', 'FontWeight', 'bold', 'Color', 'k'); 
    ylabel(ax1, 'Velocity (m/s)', 'FontWeight', 'bold', 'Color', 'k'); 
    grid on;
    xlim(ax1, [-0.07 0.07]); ylim(ax1, [-1.5 1.5]);
    
    ax2 = subplot(1,2,2);
    set(ax2, 'Color', [0.1 0.1 0.1], 'XColor', 'k', 'YColor', 'k');
    title(ax2, 'Pitch Phase Portrait', 'FontWeight', 'bold', 'Color', 'k');
    xlabel(ax2, 'Angle (rad)', 'FontWeight', 'bold', 'Color', 'k'); 
    ylabel(ax2, 'Angular Velocity (rad/s)', 'FontWeight', 'bold', 'Color', 'k'); 
    grid on;
    xlim(ax2, [-0.06 0.06]); ylim(ax2, [-1.5 1.5]);
    
    % Initialize animated lines (Fastest way to draw trajectories in MATLAB)
    line_z = animatedline(ax1, 'Color', 'c', 'LineWidth', 1.5);
    line_phi = animatedline(ax2, 'Color', 'm', 'LineWidth', 1.5);
    
    % The time tracker remains white ('w') because it sits on the dark gray axis background
    txt_tracker = text(ax1, -0.06, 1.3, 'Time: 0.00 s', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w');
    
    % 3. ANIMATION & GIF EXPORT LOOP
    gif_filename = 'Tractor_Phases.gif';
    frame_skip = 5;
    
    for k = 1:frame_skip:length(t)
        % Add points to the animated line
        addpoints(line_z, zc(k), zc_dot(k));
        addpoints(line_phi, phi(k), phi_dot(k));
        set(txt_tracker, 'String', sprintf('Time: %.2f s', t(k)));
        
        drawnow limitrate;
        
        % Capture and save frame to GIF
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

function dx = state_space_matrix(t, x, mc, Jc, kf, kr, cf, cr, lf, lr, v)
    delay = (lf + lr) / v; 
    freq = (10 / 10) * t; amp = 0.02; 
    hf = amp * sin(2*pi*freq*t); 
    freq_rate = 10 / 10;
    hf_dot = amp * 2*pi * (freq + freq_rate * t) * cos(2*pi*freq*t);
    if t >= delay
        freq_del = (10 / 10) * (t-delay);
        hr = amp * sin(2*pi*freq_del*(t-delay));
        hr_dot = amp * 2*pi * (freq_del + freq_rate * (t-delay)) * cos(2*pi*freq_del*(t-delay));
    else
        hr = 0; hr_dot = 0;
    end
    Fz = kf*hf + cf*hf_dot + kr*hr + cr*hr_dot;
    Mphi = cr*lr*hr_dot - cf*lf*hf_dot + kr*lr*hr - kf*lf*hf;
    dx = [x(2); (1/mc) * (Fz - (cf+cr)*x(2) - (kf+kr)*x(1) - (cr*lr - cf*lf)*x(4) - (kr*lr - kf*lf)*x(3)); ...
          x(4); (1/Jc) * (Mphi - (cr*lr - cf*lf)*x(2) - (kr*lr - kf*lf)*x(1) - (cf*lf^2 + cr*lr^2)*x(4) - (kf*lf^2 + kr*lr^2)*x(3))];
end