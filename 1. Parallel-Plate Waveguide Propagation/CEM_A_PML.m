c = 3e8;
f = 3.75e9;
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

Xn_phys = round(Xmax / dx);
Yn = round(Ymax / dy);

% PML
Npml = 8;
O = 2;
R = 1e-6;

Xn = Xn_phys + 2*Npml;

% χρόνος
dt_c = dx / (sqrt(2)*c);
dt = 0.99*dt_c;
Tmax = 10*T;
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


% video_filename = 'FDTD_waveguide_animation3.mp4';
% v = VideoWriter(video_filename, 'MPEG-4');
% v.FrameRate = 10;
% v.Quality = 95; 
% 
% open(v); 
% 
% fig = figure('Name', 'FDTD Recording', 'Color', 'w');
% axis tight; 
% colormap('jet');


figure('Name', 'FDTD TM Waveguide Propagation', 'Color', 'w');
title('Διάδοση κύματος @3.75GHz, Δt=0.99Δt_c με PML');
axis tight;
colormap('jet');


% PML-αρχικοποίηση
Hx_pml_L = zeros(Npml, Yn);
Hy_pml_L = zeros(Npml, Yn+1);
Ezx_pml_L = zeros(Npml, Yn+1);
Ezy_pml_L = zeros(Npml, Yn+1);

Hx_pml_R = zeros(Npml, Yn);
Hy_pml_R = zeros(Npml, Yn+1);
Ezx_pml_R = zeros(Npml, Yn+1);
Ezy_pml_R = zeros(Npml, Yn+1);

sigmaE = zeros(1, Npml);
sigmaHx = zeros(1, Npml);
sigmaHy = zeros(1, Npml);

se = - e0*c*log(R) / (2^(O+2)*dx*Npml^(O+1));
for i=1:Npml
    sigmaE(i) = se*((2*i+1)^(O+1)-(2*i-1)^(O+1));
end

sh = m0/e0 *se;
for i=1:Npml
    sigmaHx(i)=sh*( (2*i+1)^(O+1) - (2*i-1)^(O+1) );
    sigmaHy(i)=sh*( (2*(i+0.5) + 1)^(O+1) - (2*(i+0.5) - 1)^(O+1) );
end
sigmaE=fliplr(sigmaE); 
sigmaHx=fliplr(sigmaHx);
sigmaHy=fliplr(sigmaHy);

% σταθερές για PML
Ca_pml = exp(1).^ (-sigmaE.*dt/e0);
Cb_pml = (1-Ca_pml) ./ (sigmaE.*dx);
Dax_pml = exp(1).^ (-sigmaHx.*dt/m0);
Dbx_pml = (1-Dax_pml) ./ (sigmaHx.*dx);
Day_pml = exp(1).^ (-sigmaHy.*dt/m0);
Dby_pml = (1-Day_pml) ./ (sigmaHy.*dx);


% ορισμός σημείων παρατήρησης
E1 = zeros(1, Tn);
E2 = zeros(1, Tn);

x1 = Npml + round(0.7*Xn_phys);
x2 = x1 + 5;
y1 = round(0.5*Yn);
y2 = y1;


% loop ενημέρωσης πινάκων
for n = 1:Tn
    % loop Ez
    for i = (Npml+2):(Xn-Npml)
        for j = 2:Yn
            Ez(i, j) = Ca(i, j) * Ez(i, j) + Cb(i, j) * ( (Hy(i, j) - Hy(i-1, j)) - Hx(i, j) + Hx(i, j-1) );
        end
    end
    
    % loop Ez_pml αριστερά
    for i = 2:Npml
        for j = 2:Yn
            Ezx_pml_L(i, j) = Ca_pml(i-1) * Ezx_pml_L(i, j) + Cb_pml(i-1) * (Hy_pml_L(i, j) - Hy_pml_L(i-1, j));
            Ezy_pml_L(i, j) = Ca_pml(i-1) * Ezy_pml_L(i, j) + Cb_pml(i-1) * (Hx_pml_L(i, j-1) - Hx_pml_L(i, j));
            Ez(i, j) = Ezx_pml_L(i, j) + Ezy_pml_L(i, j);
        end
    end
    % άκρη του PML-κυρίως πλέγματος
    Ez(Npml+1, 2:Yn) = Ca(Npml+1, 2:Yn) .* Ez(Npml+1, 2:Yn) + Cb(Npml+1, 2:Yn) .* ( (Hy(Npml+1, 2:Yn) - Hy_pml_L(Npml, 2:Yn)) - Hx(Npml+1, 2:Yn) + Hx(Npml+1, 1:(Yn-1)) );
    % loop Ez_pml δεξιά
    for i = (Xn+2-Npml):Xn
        k = Xn + 1 - i; % αντίστροφος δείκτης για δεξί όριο
        for j = 2:Yn
            Ezx_pml_R(k, j) = Ca_pml(k) * Ezx_pml_R(k, j) + Cb_pml(k) * (Hy_pml_R(k, j) - Hy_pml_R(k+1, j));
            Ezy_pml_R(k, j) = Ca_pml(k) * Ezy_pml_R(k, j) + Cb_pml(k) * (Hx_pml_R(k, j-1) - Hx_pml_R(k, j));
            Ez(i, j) = Ezx_pml_R(k, j) + Ezy_pml_R(k, j);
        end
    end
    % άκρη κυρίως πλέγματος-PML
    Ez(Xn+1-Npml,2:Yn) = Ca(Xn+1-Npml, 2:Yn) .* Ez(Xn+1-Npml, 2:Yn) + Cb(Xn+1-Npml, 2:Yn) .* ( (Hy_pml_R(Npml, 2:Yn) - Hy(Xn-Npml, 2:Yn)) - Hx(Xn+1-Npml, 2:Yn) + Hx(Xn+1-Npml, 1:(Yn-1)) );


    % προσθήκη άμεσης πηγής στη μέση του χώρου (έστω Α0=1)
    i_source = Xn/2 + 1;
    j = 2:Yn;
    Ez(i_source, j) = 1* sin(pi*(j-1)*dy/b)*sin(2*pi*f*n*dt);


    % loop Hx
    for i = (Npml+1):(Xn-Npml+1)
        for j = 1:Yn
            Hx(i,j) = Da * Hx(i,j) - Db * (Ez(i, j+1) - Ez(i,j));
        end
    end
    % loop Hx_pml αριστερά
    for i = 2:Npml
        for j = 1:Yn
            Hx_pml_L(i, j) = Dax_pml(i-1) * Hx_pml_L(i, j) + Dbx_pml(i-1) * (Ezx_pml_L(i, j) + Ezy_pml_L(i, j) - Ezx_pml_L(i, j+1) - Ezy_pml_L(i, j+1));
        end
    end
    % loop Hx_pml δεξιά
    for i = (Xn+2-Npml):Xn
        k = Xn + 1 - i; 
        for j = 1:Yn
            Hx_pml_R(k, j) = Dax_pml(k) * Hx_pml_R(k, j) + Dbx_pml(k) * (Ezx_pml_R(k, j) + Ezy_pml_R(k, j) - Ezx_pml_R(k, j+1) - Ezy_pml_R(k, j+1));
        end
    end


    % loop Hy
    for i = (Npml+1):(Xn-Npml)
        for j = 2:Yn
            Hy(i,j) = Da * Hy(i,j) + Db * (Ez(i+1, j) - Ez(i,j));
        end
    end
    % loop Hy_pml αριστερά
    for i = 1:(Npml-1)
        for j = 2:Yn
            Hy_pml_L(i, j) = Day_pml(i) * Hy_pml_L(i, j) + Dby_pml(i) * (Ezx_pml_L(i+1, j) + Ezy_pml_L(i+1, j) - Ezx_pml_L(i, j) - Ezy_pml_L(i, j));
        end
    end
    % άκρη του pml-κυρίως πλέγματος
    Hy_pml_L(Npml, 2:Yn) = Day_pml(Npml) * Hy_pml_L(Npml, 2:Yn) + Dby_pml(Npml) * (Ez(Npml+1, 2:Yn) - Ezx_pml_L(Npml, 2:Yn) - Ezy_pml_L(Npml, 2:Yn));
    % loop Hy_pml δεξιά
    for i = (Xn+2-Npml):(Xn-1)
        k = Xn + 1 - i;
        for j = 2:Yn
            Hy_pml_R(k, j) = Day_pml(k) * Hy_pml_R(k, j) + Dby_pml(k) * (Ezx_pml_R(k-1, j) + Ezy_pml_R(k-1, j) - Ezx_pml_R(k, j) - Ezy_pml_R(k, j));
        end
    end
    % δεξιά άκρη του PML
    Hy_pml_R(1, 2:Yn) = Day_pml(1) * Hy_pml_R(1, 2:Yn) + Dby_pml(1) * (0 - Ezx_pml_R(1, 2:Yn) - Ezy_pml_R(1, 2:Yn));
    % άκρη κυρίως πλέγματος-PML
    Hy_pml_R(Npml, 2:Yn) = Day_pml(Npml) * Hy_pml_R(Npml, 2:Yn) + Dby_pml(Npml) * (Ezx_pml_R(Npml-1, 2:Yn) + Ezy_pml_R(Npml-1, 2:Yn) - Ez(Xn+1-Npml, 2:Yn));


    imagesc(Ez');
    axis equal tight;
    clim([-1, 1]);
    drawnow;

    % recording_step_interval = 3;
    % if mod(n, recording_step_interval) == 0
    % 
    %     imagesc(Ez');
    %     title('Διάδοση κύματος @3.75GHz, Δt=0.99Δt_c με PML')
    %     axis equal tight; 
    %     clim([-1, 1]); 
    %     frame = getframe(fig);
    %     writeVideo(v, frame);
    % end

    E1(n) = Ez(x1,y1);
    E2(n) = Ez(x2, y2);

end

close(v); 

disp(['Το βίντεο αποθηκεύτηκε ως: ' video_filename]);

t0 = round(7*T/dt);
t1 = round( 7.25*T /dt); % t0 + T/4

I1 = E1(t0);
Q1 = E1(t1);
I2 = E2(t0);
Q2 = E2(t1);

phi1 = atan2(Q1, I1);
phi2 = atan2(Q2, I2);


up = (2*pi*f*(x2-x1)*dx) / mod(phi2-phi1, 2*pi);
up_th = 5e8;

error = abs(up-up_th)/up_th*100;

fprintf('u_p (FDTD) = %.4e m/s\n', up);
fprintf('error = %.2f%%\n', error);