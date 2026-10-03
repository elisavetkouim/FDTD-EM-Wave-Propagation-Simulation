c = 3e8;
f = 2.5e9;
lambda = c / f;
T = 1/f;

e0 = 8.85418781762e-12;
m0 = 1.256637061436e-6;

% ζητούμενη διάσταση κ/ου
b = c / (2*3e9);

% βήμα διακριτοποίησης
dx = lambda / 32;
dy = dx;

% όρια χώρου
Xmax = 6 * lambda;
Ymax = b;

Xn = round(Xmax / dx);
Yn = round(Ymax / dy);

% χρόνος
dt_c = dx / (sqrt(2)*c);
dt = 0.99*dt_c;
Tmax = 5*T;
Tn = round(Tmax / dt);

% αρχικοποίηση πινάκων
Ez = zeros(Xn+1, Yn+1);

e(1:(Xn+1),1:(Yn+1)) = e0;
sigma = zeros(Xn+1, Yn+1);

Hx = zeros(Xn+1, Yn);
Hy = zeros(Xn, Yn+1);

m = m0;
sigma_m = 0;

% σταθερές για υπολογισμό των συνιστοσών πεδίων
Ca = ( 1 - (sigma*dt)./(2*e) ) ./ ( 1 + (sigma*dt)./(2*e) );
Cb = (dt) ./ ( e*dx.*( 1 +  (sigma*dt)./(2*e) ) );

Da = ( 1 - (sigma_m*dt)/(2*m) ) / ( 1 + (sigma_m*dt)/(2*m) );
Db = (dt) / ( m*dy*( 1 +  (sigma_m*dt)/(2*m) ) );

% video_filename = 'FDTD_waveguide_animation2c.mp4';
% v = VideoWriter(video_filename, 'MPEG-4');
% v.FrameRate = 10;
% v.Quality = 80; 
% 
% open(v); 
% 
% fig = figure('Name', 'FDTD Recording', 'Color', 'w');
% axis tight; 
% colormap('jet');

figure('Name', 'FDTD TM Waveguide Propagation', 'Color', 'w');
title('Διάδοση κύματος @2.5GHz, Δt=0.99Δt_c');
axis tight;
colormap('jet');

% loop ενημέρωσης πινάκων
for n = 1:Tn
    % loop Ez
    for i = 2:Xn
        for j = 2:Yn
            Ez(i, j) = Ca(i, j) * Ez(i, j) + Cb(i, j) * ( (Hy(i, j) - Hy(i-1, j)) - Hx(i, j) + Hx(i, j-1) );
        end
    end
    
    % προσθήκη άμεσης πηγής στη μέση του χώρου (έστω Α0=1)
    i_source = Xn/2 + 1;
    j = 2:Yn;
    Ez(i_source, j) = 1* sin(pi*(j-1)*dy/b)*sin(2*pi*f*n*dt);
    
    % loop Hx
    for i = 1:Xn
        for j = 1:Yn
            Hx(i,j) = Da * Hx(i,j) - Db * (Ez(i, j+1) - Ez(i,j));
        end
    end
    
    % loop Hy
    for i = 1:Xn
        for j = 2:Yn
            Hy(i,j) = Da * Hy(i,j) + Db * (Ez(i+1, j) - Ez(i,j));
        end
    end

    % recording_step_interval = 3;
    % if mod(n, recording_step_interval) == 0
    % 
    %     imagesc(Ez');
    %     title('Διάδοση κύματος @2.5GHz, Δt=0.99Δt_c')
    %     axis equal tight; 
    %     clim([-1, 1]); 
    %     frame = getframe(fig);
    %     writeVideo(v, frame);
    % end

    imagesc(Ez');
    axis equal tight;
    clim([-1, 1]);
    drawnow;

end
% close(v); 
% 
% disp(['Το βίντεο αποθηκεύτηκε ως: ' video_filename]);