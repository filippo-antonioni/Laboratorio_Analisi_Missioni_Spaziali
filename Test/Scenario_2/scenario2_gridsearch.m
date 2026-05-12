% =========================================================================
% SCENARIO #2: Trasferimento Diretto Terra -> Asteroide 363505 (2003 UC20)
% Metodo: Grid-Search Iterativa con Restringimento Automatico
% =========================================================================
clear; clc; close all;

% --- 1. DATI E COSTANTI FISICHE ---
mu_sun  = 1.32712440018e11; % Parametro gravitazionale del Sole [km^3/s^2]
R_sun   = 696340;           % Raggio del Sole [km]
AU      = 149597870.7;      % 1 Unità Astronomica [km]
margine = 50000000;         % Margine di sicurezza per evitare il Sole [km]

% Dati Terra
a_T  = 1.4946e8;     
e_T  = 0.016;        
i_T  = 9.1920e-5;    
OM_T = 2.7847;       
om_T = 5.2643;

% Dati Asteroide 23: 363505 (2003 UC20) 
a_A  = 0.781241 * AU;       
e_A  = 0.336932;            
i_A  = 3.78 * (pi/180);     
OM_A = 187.92 * (pi/180);   
om_A = 60.16 * (pi/180);    


% --- 2. PARAMETRI DI OTTIMIZZAZIONE ITERATIVA ---
tolleranza_DV = 0.0001; % Tolleranza per fermare la ricerca [km/s]
max_iter      = 10;     % Numero massimo di iterazioni (zoom-in) per sicurezza
N_punti       = 301;    % Punti per griglia (numero dispari per convergenza monotona)

% Centri iniziali e ampiezze (partiamo esplorando tutto il cerchio da 0 a 2*pi)
centro_th1 = pi; ampiezza_th1 = pi; 
centro_th2 = pi; ampiezza_th2 = pi;
centro_omT = pi; ampiezza_omT = pi;

diff_DV = inf;       % Inizializza la differenza con un valore enorme
DV_min_old = 1e6;    % Valore fittizio di partenza
iter = 1;            % Contatore cicli

fprintf('Inizio Ottimizzazione Iterativa...\n');
fprintf('Tolleranza impostata: %.4f km/s\n\n', tolleranza_DV);

% --- 3. CICLO WHILE (RESTRINGIMENTO AUTOMATICO) ---
storia_DV = [];

while diff_DV > tolleranza_DV && iter <= max_iter
    
    fprintf('--- Macro-Iterazione %d ---\n', iter);
    
    % Creazione dei vettori griglia limitati intorno al centro attuale
    theta1_vec = linspace(centro_th1 - ampiezza_th1, centro_th1 + ampiezza_th1, N_punti);
    theta2_vec = linspace(centro_th2 - ampiezza_th2, centro_th2 + ampiezza_th2, N_punti);
    omegaT_vec = linspace(centro_omT - ampiezza_omT, centro_omT + ampiezza_omT, N_punti);
    
    DV_min = inf; % Resetta il minimo locale per questa griglia
    
    for i = 1:length(theta1_vec)
        th1 = theta1_vec(i);
        [r1, v1i] = par2car(a_T, e_T, i_T, OM_T, om_T, th1, mu_sun);
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
                omT = omegaT_vec(k);
                
                R_om = [ cos(omT),  sin(omT), 0;
                        -sin(omT),  cos(omT), 0;
                               0,         0,  1];
                
                T_Elio_PF = R_om * R_i * R_OM; 
                
                r1_PF = T_Elio_PF * r1;
                r2_PF = T_Elio_PF * r2;
                
                th1_T = atan2(r1_PF(2), r1_PF(1));
                th2_T = atan2(r2_PF(2), r2_PF(1));
                
                num_e = norm_r2 - norm_r1;
                den_e = norm_r1 * cos(th1_T) - norm_r2 * cos(th2_T);
                
                if abs(den_e) < 1e-6
                    continue; 
                end
                
                e_T_calc = num_e / den_e;
                
                % VINCOLI FISICI
                if e_T_calc < 0 || e_T_calc >= 1
                    continue;
                end
                
                p_T_calc = norm_r1 * (1 + e_T_calc * cos(th1_T));
                a_T_calc = p_T_calc / (1 - e_T_calc^2);
                
                r_peri = a_T_calc * (1 - e_T_calc);
                if r_peri <= (R_sun + margine)
                    continue;
                end
                
                % CALCOLO DELLE VELOCITÀ tramite par2car
                [~, v1T] = par2car(a_T_calc, e_T_calc, i_trasf, OM_transf, omT, th1_T, mu_sun);
                [~, v2T] = par2car(a_T_calc, e_T_calc, i_trasf, OM_transf, omT, th2_T, mu_sun);
                
                DV_1 = norm(v1T - v1i);
                DV_2 = norm(v2f - v2T);
                DV_tot = DV_1 + DV_2;
                
                % Aggiornamento minimo locale
                if DV_tot < DV_min
                    DV_min = DV_tot;
                    ottimo_th1 = th1;
                    ottimo_th2 = th2;
                    ottimo_omT = omT;
                    ottimo_eT  = e_T_calc;
                    ottimo_aT  = a_T_calc;
                end
            end % Fine ciclo omega_T
        end % Fine ciclo theta_2
    end % Fine ciclo theta_1
    
    % --- AGGIORNAMENTO PER L'ITERAZIONE SUCCESSIVA ---
    diff_DV = abs(DV_min_old - DV_min);
    fprintf('DV trovato: %.4f km/s (Miglioramento: %.4f km/s)\n', DV_min, diff_DV);
    
    DV_min_old = DV_min;
    storia_DV = [storia_DV, DV_min];
    
    % I nuovi centri diventano gli ottimi appena trovati
    centro_th1 = ottimo_th1;
    centro_th2 = ottimo_th2;
    centro_omT = ottimo_omT;
    
    % Riduciamo l'ampiezza di ricerca del 50%
    ampiezza_th1 = ampiezza_th1 / 2;
    ampiezza_th2 = ampiezza_th2 / 2;
    ampiezza_omT = ampiezza_omT / 2;
    
    iter = iter + 1;
end

% --- 4. STAMPA DEI RISULTATI FINALI ---
if isinf(DV_min)
    fprintf('\nNessuna orbita di trasferimento fisicamente possibile trovata.\n');
else
    fprintf('\n==========================================\n');
    fprintf('--- SOLUZIONE OTTIMA GLOBALE TROVATA ---\n');
    fprintf('==========================================\n');
    fprintf('Iterazioni completate  : %d\n', iter-1);
    fprintf('Delta V Totale Minimo  : %.4f km/s\n', DV_min);
    
    % Riportiamo gli angoli tra 0 e 360 gradi usando mod()
    fprintf('Anomalia vera Terra    : %.2f deg\n', rad2deg(mod(ottimo_th1, 2*pi)));
    fprintf('Anomalia vera Asteroide: %.2f deg\n', rad2deg(mod(ottimo_th2, 2*pi)));
    fprintf('Argomento pericentro T.: %.2f deg\n', rad2deg(mod(ottimo_omT, 2*pi)));
    
    fprintf('Eccentricità trasf.    : %.4f\n', ottimo_eT);
    fprintf('Semiasse magg. trasf.  : %.2f km (%.3f AU)\n', ottimo_aT, ottimo_aT/AU);
end

% =========================================================================
%  --- 5. PLOTTING DEL PROCESSO E DEL RISULTATO ---
%  =========================================================================
if ~isinf(DV_min)
    
    % ---------------------------------------------------------
    % FIGURA 1: Convergenza del Processo (Grid-Search Zoom-in)
    % ---------------------------------------------------------
    if exist('storia_DV', 'var')
        figure('Name', 'Processo di Ottimizzazione');
        plot(1:length(storia_DV), storia_DV, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
        grid on;
        title('Convergenza del \DeltaV Totale Minimo (Grid Search)');
        xlabel('Macro-Iterazioni (Zoom-in)');
        ylabel('\DeltaV Totale [km/s]');
        xlim([1, length(storia_DV)]);
    end
    
    % ---------------------------------------------------------
    % RICALCOLO PARAMETRI ORBITA DI TRASFERIMENTO OTTIMA
    % ---------------------------------------------------------
    [r1_opt, ~] = par2car(a_T, e_T, i_T, OM_T, om_T, ottimo_th1, mu_sun);
    [r2_opt, ~] = par2car(a_A, e_A, i_A, OM_A, om_A, ottimo_th2, mu_sun);
    
    % Ricalcolo Inclinazione (i) e RAAN (OM) del piano di trasferimento
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
    
    % ---------------------------------------------------------
    % FIGURA 2: Visualizzazione delle Orbite 3D
    % ---------------------------------------------------------
    figure('Name', 'Orbite nel Sistema Solare');
    hold on; grid on; axis equal; view(3);
    
    % Il Sole a scala reale (con bordo scuro a contrasto)
    plot3(0, 0, 0, 'oy', 'MarkerSize', 8, 'MarkerFaceColor', '#FFCC00', 'MarkerEdgeColor', 'y', 'LineWidth', 1.5, 'DisplayName', 'Sole');
    
    % Generazione punti per tracciare le ellissi complete
    theta_plot = linspace(0, 2*pi, 300);
    r_T_plot = zeros(3, 300);
    r_A_plot = zeros(3, 300);
    r_Trasf_plot = zeros(3, 300);
    
    for idx = 1:length(theta_plot)
        [r_T_plot(:,idx), ~] = par2car(a_T, e_T, i_T, OM_T, om_T, theta_plot(idx), mu_sun);
        [r_A_plot(:,idx), ~] = par2car(a_A, e_A, i_A, OM_A, om_A, theta_plot(idx), mu_sun);
        [r_Trasf_plot(:,idx), ~] = par2car(ottimo_aT, ottimo_eT, i_trasf_opt, OM_transf_opt, ottimo_omT, theta_plot(idx), mu_sun);
    end
    
    % Plot Orbita Terra
    plot3(r_T_plot(1,:), r_T_plot(2,:), r_T_plot(3,:), 'b', 'LineWidth', 1.2, 'DisplayName', 'Orbita Terra');
    
    % Plot Orbita Asteroide
    plot3(r_A_plot(1,:), r_A_plot(2,:), r_A_plot(3,:), 'r', 'LineWidth', 1.2, 'DisplayName', 'Orbita Asteroide 363505');
    
    % Plot Orbita di Trasferimento (intera, tratteggiata)
    plot3(r_Trasf_plot(1,:), r_Trasf_plot(2,:), r_Trasf_plot(3,:), '--g', 'LineWidth', 1.5, 'DisplayName', 'Orbita Trasferimento');
    
    % Plot Punti di Partenza (Terra) e Arrivo (Asteroide)
    plot3(r1_opt(1), r1_opt(2), r1_opt(3), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', 'Partenza (Terra)');
    plot3(r2_opt(1), r2_opt(2), r2_opt(3), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Arrivo (Asteroide)');
    
    % --- EVIDENZIAZIONE DEL TRATTO EFFETTIVAMENTE PERCORSO ---
    R_OM_opt = [ cos(OM_transf_opt),  sin(OM_transf_opt), 0;
                -sin(OM_transf_opt),  cos(OM_transf_opt), 0;
                       0,               0,        1];
    R_i_opt =  [1,        0,               0;
                0,  cos(i_trasf_opt),  sin(i_trasf_opt);
                0, -sin(i_trasf_opt),  cos(i_trasf_opt)];
    R_om_opt = [ cos(ottimo_omT),  sin(ottimo_omT), 0;
                -sin(ottimo_omT),  cos(ottimo_omT), 0;
                       0,         0,  1];
    T_Elio_PF_opt = R_om_opt * R_i_opt * R_OM_opt; 
    
    r1_PF_opt = T_Elio_PF_opt * r1_opt;
    r2_PF_opt = T_Elio_PF_opt * r2_opt;
    
    th1_T_opt = mod(atan2(r1_PF_opt(2), r1_PF_opt(1)), 2*pi);
    th2_T_opt = mod(atan2(r2_PF_opt(2), r2_PF_opt(1)), 2*pi);
    
    if th2_T_opt < th1_T_opt
        th2_T_opt = th2_T_opt + 2*pi;
    end
    
    theta_arco = linspace(th1_T_opt, th2_T_opt, 150);
    r_Arco_plot = zeros(3, length(theta_arco));
    
    for idx = 1:length(theta_arco)
        [r_Arco_plot(:,idx), ~] = par2car(ottimo_aT, ottimo_eT, i_trasf_opt, OM_transf_opt, ottimo_omT, theta_arco(idx), mu_sun);
    end
    
    % Plottiamo l'arco sopra la linea tratteggiata (spessore 3.5, colore verde continuo)
    plot3(r_Arco_plot(1,:), r_Arco_plot(2,:), r_Arco_plot(3,:), '-g', 'LineWidth', 3.5, 'DisplayName', 'Tratto Percorso (Volo)');
    % ---------------------------------------------------------
    
    % Formattazione grafico finale
    title(sprintf('Trasferimento Diretto: Terra -> Asteroide 363505\n\\DeltaV = %.4f km/s', DV_min));
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    legend('Location', 'best');
end