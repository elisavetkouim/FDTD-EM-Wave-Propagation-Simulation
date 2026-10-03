c = 3e8;
f0 = 1.5e9;
T = 1 / f0;
lambda = c / f0;
L = 1;
Ly = 0.1;
er = 5.25;
de = 0.85;
omega = 2*pi*1e9;

e0 = 8.85418781762e-12;
m0 = 1.256637061436e-6;

% βήμα διακριτοποίησης κατα y
dy = lambda / 16;
Yn = round(2*L / dy);

% βήμα διακριτοποίησης χρόνου
dt_c = dy / c;
dt = 0.99*dt_c;
Tmax = 10*T;
Tn = round(Tmax / dt);

% αρχικοποίηση πινάκων πεδίου
Ez = zeros(1, Yn+1);
Hx = zeros(1, Yn);


% διηλεκτρική σταθερά

e(1:(Tn+1), 1:(Yn+1)) = e0; % (Tn+1) x (Yn+1)
% στις θέσεις που είναι η πλάκα
i_start = round(Yn/2+1);
i_end = i_start + round(Ly/dy);

er_t = er + de*sin(omega*dt*(0:Tn)); % 1 x (Tn+1)
er_t = er_t(:); % (Tn+1) x 1
e(:, i_start:i_end) = repmat(e0 * er_t, 1, i_end-i_start+1);

% σταθερές για υπολογισμό των συνιστοσών πεδίων
Ca = e(1:Tn, 2:Yn) ./ e(2:(Tn+1), 2:Yn);
Cb = dt ./ (dy * e(2:(Tn+1), 2:Yn) );


% video_filename = 'fdtd_simulation_no_PML.mp4';
% v = VideoWriter(video_filename, 'MPEG-4');
% v.FrameRate = 20;
% v.Quality = 95;
% open(v);
% 
% fig = figure('Name', 'FDTD Recording', 'Color', 'w');
% axis tight; 
% colormap('jet');

figure('Color', 'w');
colormap('jet');

% loop ενημέρωσης πινάκων
for n = 1:Tn

    Ez(2:Yn) = Ca(n, :) .* Ez(2:Yn) - Cb(n, :) .* ( Hx(2:Yn) - Hx(1:(Yn-1)) );
    
    % source
    i_source = 1 + round(2*lambda/dy);
    n_source_end = round(3*T/dt);
    if n <= n_source_end
        Ez(i_source) = sin(2*pi*f0*n*dt);
    end

    Hx(1:Yn) = Hx(1:Yn) - (dt/(m0*dy)) .* ( Ez(2:(Yn+1)) - Ez(1:Yn) );

    


    imagesc(Ez);
    clim([-1, 1]);
    title('Διάδοση κύματος χωρίς PML')
    axis normal;
    hold on;

    patch('XData', [i_start, i_end, i_end, i_start], ...
          'YData', [0.5, 0.5, 1.5, 1.5], ...
          'FaceColor', [1, 0, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    xline(i_source, 'g-', 'LineWidth', 2);

    drawnow;

    % recording_step_interval = 1;
    % if mod(n, recording_step_interval) == 0
    % 
    %     imagesc(Ez);
    %     clim([-1, 1]);
    %     title('Διάδοση κύματος χωρίς PML')
    %     axis normal;
    %     hold on;
    % 
    %     patch('XData', [i_start, i_end, i_end, i_start], ...
    %           'YData', [0.5, 0.5, 1.5, 1.5], ...
    %           'FaceColor', [1, 0, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    % 
    %     xline(i_source, 'g-', 'LineWidth', 2);
    %         frame = getframe(fig);
    %     writeVideo(v, frame);
    % end

    
    
end

% close(v);
% disp(['Το βίντεο αποθηκεύτηκε ως: ' video_filename]);