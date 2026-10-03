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
Yn_physical = round(2*L / dy);

% βήμα διακριτοποίησης χρόνου
dt_c = dy / c;
dt = 0.99*dt_c;
Tmax = 16*T;
Tn = round(Tmax / dt);

% PML
Npml = 100;
O = 2;
R = 1e-7;

Yn = Yn_physical + 2*Npml;

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

% σταθερές για υπολογισμό των συνιστοσών πεδίων (μέσο χωρίς απώλειες -> σ,σ* = 0)
Ca = e(1:Tn, 2:Yn) ./ e(2:(Tn+1), 2:Yn);
Cb = dt ./ (dy * e(2:(Tn+1), 2:Yn) );


% PML-αρχικοποίηση
E_pml_L = zeros(1, Npml);
E_pml_R = zeros(1, Npml);
H_pml_L = zeros(1, Npml);
H_pml_R = zeros(1, Npml);



se = - e0*c*log(R) / (2^(O+2)*dy*Npml^(O+1));
sigmaE = se* ( (2*(1:Npml) + 1).^(O+1) -(2*(1:Npml) - 1).^(O+1) );

sh = m0*se/e0;
sigmaH = sh* ( (2*((1:Npml)+0.5) + 1).^(O+1) - (2*((1:Npml)+0.5) - 1).^(O+1) );

sigmaE = fliplr(sigmaE);
sigmaH = fliplr(sigmaH);

% σταθερές για PML
Ca_pml = exp(1).^ (-sigmaE.*dt/e0);
Cb_pml = (1 - Ca_pml) ./ (sigmaE.*dy);
Da_pml = exp(1).^ (-sigmaH.*dt/m0);
Db_pml = (1 - Da_pml) ./ (sigmaH.*dy);

Ca_pml_R = fliplr(Ca_pml);
Cb_pml_R = fliplr(Cb_pml);
Da_pml_R = fliplr(Da_pml);
Db_pml_R = fliplr(Db_pml);



video_filename = 'fdtd_simulation_PML_100.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 10;
v.Quality = 95;
open(v);

fig = figure('Name', 'FDTD Recording', 'Color', 'w');
axis tight; 
colormap('jet');

% figure('Color', 'w');
% colormap('jet');




% loop ενημέρωσης πινάκων
for n = 1:Tn
    % Ez
    Ez((Npml+2):(Yn-Npml)) = Ca(n, (Npml+2):(Yn-Npml)) .* Ez((Npml+2):(Yn-Npml)) - Cb(n, (Npml+2):(Yn-Npml)) .* ( Hx((Npml+2):(Yn-Npml)) - Hx((Npml+1):(Yn-Npml-1)) );
    % αριστερό PML
    E_pml_L(2:Npml) = Ca_pml(1:(Npml-1)) .* E_pml_L(2:Npml) - Cb_pml(1:(Npml-1)) .* (H_pml_L(2:Npml) - H_pml_L(1:(Npml-1)));
    Ez(2:Npml) = E_pml_L(2:Npml);
    % άκρη PML-κύριου πλέγματος
    Ez(Npml+1) = Ca(n, Npml+1) * Ez(Npml+1) - Cb(n, Npml+1) * (Hx(Npml+1) - H_pml_L(Npml));
    % δεξί PML
    E_pml_R(1:(Npml-1)) = Ca_pml_R(1:(Npml-1)).* E_pml_R(1:(Npml-1)) - Cb_pml_R(1:(Npml-1)) .* (H_pml_R(2:Npml) - H_pml_R(1:(Npml-1)));
    Ez((Yn+2-Npml):Yn) = E_pml_R(1:(Npml-1));
    % κύριο πλέγμα-άκρη PML
    Ez(Yn+1-Npml) = Ca(n, Yn+1-Npml) * Ez(Yn+1-Npml) - Cb(n, Yn+1-Npml) * (H_pml_R(1) - Hx(Yn-Npml));
    
    % source
    i_source = Npml + 1 + round(2*lambda/dy);
    n_source_end = round(3*T/dt);
    if n <= n_source_end
        Ez(i_source) = sin(2*pi*f0*n*dt);
    end

    % Hx
    Hx((Npml+1):(Yn-Npml)) = Hx((Npml+1):(Yn-Npml)) - (dt/(m0*dy)) .* ( Ez((Npml+2):(Yn-Npml+1)) - Ez((Npml+1):(Yn-Npml)) );
    % αριστερό PML
    H_pml_L(1:(Npml-1)) = Da_pml(1:(Npml-1)) .* H_pml_L(1:(Npml-1)) - Db_pml(1:(Npml-1)) .* (E_pml_L(2:(Npml)) - E_pml_L(1:(Npml-1)));
    % άκρη PML-κύριου πλέγματος
    H_pml_L(Npml) = Da_pml(Npml) * H_pml_L(Npml) - Db_pml(Npml) * (Ez(Npml+1) - E_pml_L(Npml));
    % δεξί PML
    H_pml_R(2:Npml) = Da_pml_R(2:Npml) .* H_pml_R(2:Npml) - Db_pml_R(2:Npml) .* (E_pml_R(2:Npml) - E_pml_R(1:(Npml-1)));
    % κύριο πλέγμα-δεξί PML
    H_pml_R(1) = Da_pml_R(1) * H_pml_R(1) - Db_pml_R(1) * (E_pml_R(1) - Ez(Yn+1-Npml));



    frame_stride = 2;
    if mod(n, frame_stride) == 0

        imagesc(Ez);
        clim([-1, 1]);
        title(['Διάδοση κύματος με PML ', num2str(Npml), ' κελιών']);
        axis normal;
        colorbar;
        hold on;
        
        patch('XData', [i_start, i_end, i_end, i_start], ...
              'YData', [0.5, 0.5, 1.5, 1.5], ...
              'FaceColor', [1, 0, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        
        xline(Npml + 1 + round(2*lambda/dy), 'g-', 'LineWidth', 2);
        hold off;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    
   
    

    % imagesc(Ez);
    % clim([-1, 1]);
    % title('Διάδοση κύματος χωρίς PML')
    % axis normal;
    % hold on;
    % 
    % patch('XData', [i_start, i_end, i_end, i_start], ...
    %       'YData', [0.5, 0.5, 1.5, 1.5], ...
    %       'FaceColor', [1, 0, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    % 
    % xline(i_source, 'g-', 'LineWidth', 2);
    % 
    % drawnow;

    
end

close(v);
disp(['Το βίντεο αποθηκεύτηκε ως: ' video_filename]);