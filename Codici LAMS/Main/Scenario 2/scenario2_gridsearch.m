clear; 
clc;
close all;

mu_sun  = 1.32712440018e11; % Parametro gravitazionale del Sole [km^3/s^2]
R_sun   = 696340;           % Raggio del Sole [km]
AU      = 149597870.7;      % 1 Unità Astronomica [km]
margine = 50000000;         % Margine di sicurezza per evitare il Sole [km]

% Dati Terra
a_E  = 1.4946e8;     
e_T  = 0.016;        
i_E  = 9.1920e-5;    
OM_E = 2.7847;       
om_E = 5.2643;

% Dati Asteroide
a_A  = 0.781241 * AU;       
e_A  = 0.336932;            
i_A  = 3.78 * (pi/180);     
OM_A = 187.92 * (pi/180);   
om_A = 60.16 * (pi/180);    


% parametri di ottimizzazione
tolleranza_DV = 0.0001; % Tolleranza per fermare la ricerca [km/s]
max_iter      = 10;     % Numero massimo di iterazioni (zoom-in) per sicurezza
N_punti       = 301;    % Punti per griglia (numero dispari per convergenza monotona)


% Centri iniziali e ampiezze (dimensione iniziale da 0 a 2pi)
centro_th1 = pi; ampiezza_th1 = pi; 
centro_th2 = pi; ampiezza_th2 = pi;
centro_om_trasf = pi; ampiezza_om_trasf = pi;

diff_DV = inf;     
DV_min_old = 1e6;  
iter = 0;       

% ciclo while con restringimento dell'intervallo
storia_DV = [];

while diff_DV > tolleranza_DV && iter <= max_iter
    iter=iter+1;
    
    theta1_vec = linspace(centro_th1 - ampiezza_th1, centro_th1 + ampiezza_th1, N_punti);
    theta2_vec = linspace(centro_th2 - ampiezza_th2, centro_th2 + ampiezza_th2, N_punti);
    omegaT_vec = linspace(centro_om_trasf - ampiezza_om_trasf, centro_om_trasf + ampiezza_om_trasf, N_punti);
    
    DV_min = inf;
    
    for i = 1:length(theta1_vec)
        th1 = theta1_vec(i);
        [r1, v1i] = par2car(a_E, e_T, i_E, OM_E, om_E, th1, mu_sun);
        norm_r1 = norm(r1);
        
        for j = 1:length(theta2_vec)
            th2 = theta2_vec(j);
            [r2, v2f] = par2car(a_A, e_A, i_A, OM_A, om_A, th2, mu_sun);
            norm_r2 = norm(r2);
            
            h_vec = cross(r1, r2);
            norm_h = norm(h_vec);
            
            if norm_h < 1e-6
                continue; 
            end
            
            h_vers = h_vec / norm_h;
            i_trasf = acos(h_vers(3)); 
            
            K = [0; 0; 1];
            N_vec = cross(K, h_vers);
            norm_N = norm(N_vec);
            
            if norm_N < 1e-6
                OM_transf = 0; 
            else
                N_vers = N_vec / norm_N;
                if N_vers(2) >= 0
                    OM_transf = acos(N_vers(1));
                else
                    OM_transf = 2*pi - acos(N_vers(1));
                end
            end
            
            R_OM = [ cos(OM_transf),  sin(OM_transf), 0;
                    -sin(OM_transf),  cos(OM_transf), 0;
                           0,               0,        1];
                   
            R_i =  [1,        0,               0;
                    0,  cos(i_trasf),  sin(i_trasf);
                    0, -sin(i_trasf),  cos(i_trasf)];
            
            for k = 1:length(omegaT_vec)
                om_trasf = omegaT_vec(k);
                
                R_om = [ cos(om_trasf),  sin(om_trasf), 0;
                        -sin(om_trasf),  cos(om_trasf), 0;
                               0,         0,  1];
                
                T_Elio_PF = R_om * R_i * R_OM; 
                
                r1_PF = T_Elio_PF * r1;
                r2_PF = T_Elio_PF * r2;
                
                th1_trasf = atan2(r1_PF(2), r1_PF(1));
                th2_trsf = atan2(r2_PF(2), r2_PF(1));
                
                num_e = norm_r2 - norm_r1;
                den_e = norm_r1 * cos(th1_trasf) - norm_r2 * cos(th2_trsf);
                
                if abs(den_e) < 1e-6
                    continue; 
                end
                
                e_trasf_calc = num_e / den_e;
                
                % VINCOLI FISICI
                if e_trasf_calc < 0 || e_trasf_calc >= 1
                    continue;
                end
                
                p_trasf_calc = norm_r1 * (1 + e_trasf_calc * cos(th1_trasf));
                a_trsf_calc = p_trasf_calc / (1 - e_trasf_calc^2);
                
                r_peri = a_trsf_calc * (1 - e_trasf_calc);
                if r_peri <= (R_sun + margine)
                    continue;
                end
                
                % CALCOLO VELOCITÀ 
                [~, v1_trasf] = par2car(a_trsf_calc, e_trasf_calc, i_trasf, OM_transf, om_trasf, th1_trasf, mu_sun);
                [~, v2_trasf] = par2car(a_trsf_calc, e_trasf_calc, i_trasf, OM_transf, om_trasf, th2_trsf, mu_sun);
                
                DV_1 = norm(v1_trasf - v1i);
                DV_2 = norm(v2f - v2_trasf);
                DV_tot = DV_1 + DV_2;
                
                % Aggiornamento minimo locale
                if DV_tot < DV_min
                    DV_min = DV_tot;
                    ottimo_th1 = th1;
                    ottimo_th2 = th2;
                    ottimo_om_trasf = om_trasf;
                    ottimo_e_trasf  = e_trasf_calc;
                    ottimo_a_trasf  = a_trsf_calc;
                end
            end % Fine ciclo omega_T
        end % Fine ciclo theta_2
    end % Fine ciclo theta_1
    
    % AGGIORNAMENTO PER ITERAZIONE SUCCESSIVA
    diff_DV = abs(DV_min_old - DV_min);
    fprintf('DV trovato: %.4f km/s (Miglioramento: %.4f km/s)\n', DV_min, diff_DV);
    
    DV_min_old = DV_min;
    storia_DV = [storia_DV, DV_min];
    
    % Aggiornamento dei centri con gli ottimi trovati 
    centro_th1 = ottimo_th1;
    centro_th2 = ottimo_th2;
    centro_om_trasf = ottimo_om_trasf;
    
    % Dimezzazione dell'ampiezza dell'intervallo
    ampiezza_th1 = ampiezza_th1 / 2;
    ampiezza_th2 = ampiezza_th2 / 2;
    ampiezza_om_trasf = ampiezza_om_trasf / 2;
    
end

% Risultati finali
if isinf(DV_min)
    fprintf('\nNessuna orbita di trasferimento fisicamente possibile trovata.\n');
else
    fprintf('\n==========================================\n');
    fprintf('--- SOLUZIONE OTTIMA GLOBALE TROVATA ---\n');
    fprintf('==========================================\n');
    fprintf('Iterazioni completate  : %d\n', iter-1);
    fprintf('Delta V Totale Minimo  : %.4f km/s\n', DV_min);
    fprintf('Anomalia vera Terra    : %.2f deg\n', rad2deg(mod(ottimo_th1, 2*pi)));
    fprintf('Anomalia vera Asteroide: %.2f deg\n', rad2deg(mod(ottimo_th2, 2*pi)));
    fprintf('Argomento pericentro T.: %.2f deg\n', rad2deg(mod(ottimo_om_trasf, 2*pi)));
    fprintf('Eccentricità trasf.    : %.4f\n', ottimo_e_trasf);
    fprintf('Semiasse magg. trasf.  : %.2f km (%.3f AU)\n', ottimo_a_trasf, ottimo_a_trasf/AU);
end

% Plot
if ~isinf(DV_min)
    
    % Convergenza
    if exist('storia_DV', 'var')
        figure('Name', 'Processo di Ottimizzazione');
        plot(1:length(storia_DV), storia_DV, '-oc', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
        grid on;
        title('Convergenza del \DeltaV Totale Minimo (Grid Search)');
        xlabel('Macro-Iterazioni (Zoom-in)');
        ylabel('\DeltaV Totale [km/s]');
        xlim([1, length(storia_DV)]);
    end
    
    [r1_opt, ~] = par2car(a_E, e_T, i_E, OM_E, om_E, ottimo_th1, mu_sun);
    [r2_opt, ~] = par2car(a_A, e_A, i_A, OM_A, om_A, ottimo_th2, mu_sun);
    
    h_vec_opt = cross(r1_opt, r2_opt);
    h_vers_opt = h_vec_opt / norm(h_vec_opt);
    i_trasf_opt = acos(h_vers_opt(3)); 
    
    N_vec_opt = cross([0; 0; 1], h_vers_opt);
    if norm(N_vec_opt) < 1e-6
        OM_transf_opt = 0;
    else
        N_vers_opt = N_vec_opt / norm(N_vec_opt);
        if N_vers_opt(2) >= 0
            OM_transf_opt = acos(N_vers_opt(1));
        else
            OM_transf_opt = 2*pi - acos(N_vers_opt(1));
        end
    end

    % Plot orbite
    figure('Name', 'Orbite nel Sistema Solare');
    hold on; grid on; axis equal; view(3);
    
    % Sole
    plot3(0, 0, 0, 'oy', 'MarkerSize', 8, 'MarkerFaceColor', '#FFCC00', 'MarkerEdgeColor', 'y', 'LineWidth', 1.5, 'DisplayName', 'Sole');
    
    % Generazione punti per tracciare le ellissi complete
    theta_plot = linspace(0, 2*pi, 300);
    r_E_plot = zeros(3, 300);
    r_A_plot = zeros(3, 300);
    r_Trasf_plot = zeros(3, 300);
    
    for idx = 1:length(theta_plot)
        [r_E_plot(:,idx), ~] = par2car(a_E, e_T, i_E, OM_E, om_E, theta_plot(idx), mu_sun);
        [r_A_plot(:,idx), ~] = par2car(a_A, e_A, i_A, OM_A, om_A, theta_plot(idx), mu_sun);
        [r_Trasf_plot(:,idx), ~] = par2car(ottimo_a_trasf, ottimo_e_trasf, i_trasf_opt, OM_transf_opt, ottimo_om_trasf, theta_plot(idx), mu_sun);
    end
    
    % Plot Orbita Terra
    plot3(r_E_plot(1,:), r_E_plot(2,:), r_E_plot(3,:), 'b', 'LineWidth', 1.2, 'DisplayName', 'Orbita Terra');
    
    % Plot Orbita Asteroide
    plot3(r_A_plot(1,:), r_A_plot(2,:), r_A_plot(3,:), 'r', 'LineWidth', 1.2, 'DisplayName', 'Orbita Asteroide 363505');
    
    % Plot Orbita di Trasferimento 
    plot3(r_Trasf_plot(1,:), r_Trasf_plot(2,:), r_Trasf_plot(3,:), '--g', 'LineWidth', 1.5, 'DisplayName', 'Orbita Trasferimento');
    
    % Plot Punti di Partenza (Terra) e Arrivo (Asteroide)
    plot3(r1_opt(1), r1_opt(2), r1_opt(3), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', 'Partenza (Terra)');
    plot3(r2_opt(1), r2_opt(2), r2_opt(3), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Arrivo (Asteroide)');
    
    R_OM_opt = [ cos(OM_transf_opt),  sin(OM_transf_opt), 0;
                -sin(OM_transf_opt),  cos(OM_transf_opt), 0;
                       0,               0,        1];
    R_i_opt =  [1,        0,               0;
                0,  cos(i_trasf_opt),  sin(i_trasf_opt);
                0, -sin(i_trasf_opt),  cos(i_trasf_opt)];
    R_om_opt = [ cos(ottimo_om_trasf),  sin(ottimo_om_trasf), 0;
                -sin(ottimo_om_trasf),  cos(ottimo_om_trasf), 0;
                       0,         0,  1];
    T_Elio_PF_opt = R_om_opt * R_i_opt * R_OM_opt; 
    
    r1_PF_opt = T_Elio_PF_opt * r1_opt;
    r2_PF_opt = T_Elio_PF_opt * r2_opt;
    
    th1_trasf_opt = mod(atan2(r1_PF_opt(2), r1_PF_opt(1)), 2*pi);
    th2_trasf_opt = mod(atan2(r2_PF_opt(2), r2_PF_opt(1)), 2*pi);
    
    if th2_trasf_opt < th1_trasf_opt
        th2_trasf_opt = th2_trasf_opt + 2*pi;
    end
    
    theta_arco = linspace(th1_trasf_opt, th2_trasf_opt, 150);
    r_Arco_plot = zeros(3, length(theta_arco));
    
    for idx = 1:length(theta_arco)
        [r_Arco_plot(:,idx), ~] = par2car(ottimo_a_trasf, ottimo_e_trasf, i_trasf_opt, OM_transf_opt, ottimo_om_trasf, theta_arco(idx), mu_sun);
    end
    
    plot3(r_Arco_plot(1,:), r_Arco_plot(2,:), r_Arco_plot(3,:), '-g', 'LineWidth', 3.5, 'DisplayName', 'Tratto Percorso (Volo)');
    
    title(sprintf('Trasferimento Diretto: Terra -> Asteroide 363505\n\\DeltaV = %.4f km/s', DV_min));
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    legend('Location', 'best');

    % Figura: fascio di orbite
    figure('Name', 'Fascio di Orbite Secanti (\omega_T variabile)');
    hold on; grid on; axis equal; view(3);
    
    [xS, yS, zS] = sphere(100);
    R_sun_plot = R_sun * 4;
    surf(xS * R_sun_plot, yS * R_sun_plot, zS * R_sun_plot, ...
        'FaceColor', '#FF8C00', 'EdgeColor', 'none', ...
        'FaceLighting', 'none', 'AmbientStrength', 1, 'HandleVisibility', 'off');
    plot3(0, 0, 0, 'oy', 'MarkerSize', 8, 'MarkerFaceColor', '#FF8C00', 'LineStyle', 'none', 'DisplayName', 'Sole');
    plot3(r1_opt(1), r1_opt(2), r1_opt(3), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'DisplayName', 'Punto Partenza');
    plot3(r2_opt(1), r2_opt(2), r2_opt(3), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'Punto Arrivo');
    
    step_w = deg2rad(2);
    omT_range = 0 : step_w : 2*pi;
    
    colore_azzurrino = [0.6, 0.8, 1];
    first_dashed_plotted = false;
    
    for idx_w = 1:length(omT_range)
        omT_test = omT_range(idx_w);
        
        R_OM = [ cos(OM_transf_opt),  sin(OM_transf_opt), 0;
                -sin(OM_transf_opt),  cos(OM_transf_opt), 0;
                       0,               0,        1];
        R_i =  [1,        0,               0;
                0,  cos(i_trasf_opt),  sin(i_trasf_opt);
                0, -sin(i_trasf_opt),  cos(i_trasf_opt)];
        R_om = [ cos(omT_test),  sin(omT_test), 0;
                -sin(omT_test),  cos(omT_test), 0;
                       0,         0,  1];
        
        T_Elio_PF = R_om * R_i * R_OM; 
        r1_PF = T_Elio_PF * r1_opt;
        r2_PF = T_Elio_PF * r2_opt;
        
        th1_trasf = atan2(r1_PF(2), r1_PF(1));
        th2_trsf = atan2(r2_PF(2), r2_PF(1));
        
        norm_r1 = norm(r1_opt);
        norm_r2 = norm(r2_opt);
        
        den_e = norm_r1 * cos(th1_trasf) - norm_r2 * cos(th2_trsf);
        if abs(den_e) < 1e-6
            continue; 
        end
        
        e_T_test = (norm_r2 - norm_r1) / den_e;
        
        % vincoli fisici
        if e_T_test >= 0 && e_T_test < 1
            p_T_test = norm_r1 * (1 + e_T_test * cos(th1_trasf));
            a_T_test = p_T_test / (1 - e_T_test^2);
            r_peri = a_T_test * (1 - e_T_test);
            
            if r_peri > (R_sun + margine)
                
                % Calcolo dei punti dell'orbita grigia
                theta_plot = linspace(0, 2*pi, 150);
                r_bundle = zeros(3, length(theta_plot));
                for p_idx = 1:length(theta_plot)
                    [r_bundle(:,p_idx), ~] = par2car(a_T_test, e_T_test, i_trasf_opt, OM_transf_opt, omT_test, theta_plot(p_idx), mu_sun);
                end
                
               
                if ~first_dashed_plotted
                    plot3(r_bundle(1,:), r_bundle(2,:), r_bundle(3,:), '--', 'Color', colore_azzurrino, 'LineWidth', 0.5, 'DisplayName', 'Fascio sub-ottimale (\Delta\omega_T = 2°)');
                    first_dashed_plotted = true;
                else
                    plot3(r_bundle(1,:), r_bundle(2,:), r_bundle(3,:), '--', 'Color', colore_azzurrino, 'LineWidth', 0.5, 'HandleVisibility', 'off');
                end
            end
        end
    end
    
    theta_plot = linspace(0, 2*pi, 300);
    r_Trasf_opt_plot = zeros(3, length(theta_plot));
    for p_idx = 1:length(theta_plot)
        [r_Trasf_opt_plot(:,p_idx), ~] = par2car(ottimo_a_trasf, ottimo_e_trasf, i_trasf_opt, OM_transf_opt, ottimo_om_trasf, theta_plot(p_idx), mu_sun);
    end
    
    plot3(r_Trasf_opt_plot(1,:), r_Trasf_opt_plot(2,:), r_Trasf_opt_plot(3,:), '-r', 'LineWidth', 3, 'DisplayName', 'Trasferimento OTTIMO');
    title('Fascio di Orbite fissati $\theta_i$ e $\theta_f$', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    legend('Location', 'northeast');
    hold off;
end


% TOF
TOF_sec = TOF(ottimo_a_trasf, ottimo_e_trasf, th1_trasf_opt, th2_trasf_opt, mu_sun);
TOF_giorni = TOF_sec / (24 * 3600);