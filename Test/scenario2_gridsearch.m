% =========================================================================
% SCENARIO #2: Trasferimento Diretto Terra -> Asteroide 363505 (2003 UC20)
% Metodo: Grid-Search Iterativa con Restringimento Automatico
% =========================================================================
clear; clc; close all;

% --- 1. DATI E COSTANTI FISICHE ---
mu_sun  = 1.32712440018e11; % Parametro gravitazionale del Sole [km^3/s^2]
R_sun   = 696340;           % Raggio del Sole [km]
AU      = 149597870.7;      % 1 Unità Astronomica [km]
margine = 50000000;           % Margine di sicurezza per evitare il Sole [km] -> Per evitare sublimazione dei materiali

% Dati Terra (da Scenario #2 slide 5)
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
max_iter      = 10;    % Numero massimo di iterazioni (zoom-in) per sicurezza
N_punti       = 300;    % Punti per griglia (300^3 = ~27.000.000 iterazioni a ciclo)

% Centri iniziali e ampiezze (partiamo esplorando tutto il cerchio da 0 a 2*pi)
centro_th1 = pi; ampiezza_th1 = pi; 
centro_th2 = pi; ampiezza_th2 = pi;
centro_omT = pi; ampiezza_omT = pi;

diff_DV = inf;       % Inizializza la differenza con un valore enorme
DV_min_old = 1e6;    % Valore fittizio di partenza per il DV del ciclo precedente
iter = 1;            % Contatore cicli

fprintf('Inizio Ottimizzazione Iterativa...\n');
fprintf('Tolleranza impostata: %.4f km/s\n\n', tolleranza_DV);

% --- 3. CICLO WHILE (RESTRINGIMENTO AUTOMATICO) ---
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
                
                % CALCOLO DELLE VELOCITÀ
                % Questo va verificato
                v1T_PF = sqrt(mu_sun / p_T_calc) * [-sin(th1_T); e_T_calc + cos(th1_T); 0];
                v2T_PF = sqrt(mu_sun / p_T_calc) * [-sin(th2_T); e_T_calc + cos(th2_T); 0];
                
                % Questo va verificato
                v1T = T_Elio_PF' * v1T_PF;
                v2T = T_Elio_PF' * v2T_PF;
                
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
    
    % I nuovi centri diventano gli ottimi appena trovati
    centro_th1 = ottimo_th1;
    centro_th2 = ottimo_th2;
    centro_omT = ottimo_omT;
    
    % Riduciamo l'ampiezza di ricerca (es. restringiamo del 50% il campo visivo)
    % N.B.: Puoi modificare questo fattore (es. /3 o /4) per stringere più o meno velocemente
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