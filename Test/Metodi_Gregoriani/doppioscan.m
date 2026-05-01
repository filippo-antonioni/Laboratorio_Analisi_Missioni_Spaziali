% =========================================================================
%   TRASFERIMENTO ORBITALE - SEQUENZA 5 MANOVRE 
%   OPZIONE 2: ORBITA AUSILIARIA ELLITTICA (Ottimizzazione 2D)
% =========================================================================
clear; 
clc; 
close all;

% -------------------------------------------------------------------------
% 1. CARICAMENTO DATI
% -------------------------------------------------------------------------
try
    T  = load("DatiSC1-2026.txt");
    gruppo = 23;
    
    % --- Orbita target / finale ---
    rr_f   = T(gruppo, 8:10)';
    vv_f   = T(gruppo, 11:13)';
    mu = 398600;  % [km^3/s^2]
    [a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);
    
    % --- Orbita iniziale ---
    orbita_iv = T(gruppo,2:7);
    a_i  = orbita_iv(1);
    e_i  = orbita_iv(2);
    i_i  = orbita_iv(3);
    OM_i = orbita_iv(4);
    om_i = orbita_iv(5);
    th_i = orbita_iv(6);
catch
    warning('File dati non trovato. Utilizzo dati di default per il test.');
    mu = 398600;
    a_i = 8000; e_i = 0.1; i_i = deg2rad(20); OM_i = deg2rad(30); om_i = deg2rad(40); th_i = 0;
    a_f = 25000; e_f = 0.3; i_f = deg2rad(60); OM_f = deg2rad(90); om_f = deg2rad(120); th_f = pi;
end

% -------------------------------------------------------------------------
% 2. OTTIMIZZAZIONE 2D ORBITA AUSILIARIA (a_aux e r_p_aux)
% -------------------------------------------------------------------------
r_p_i = a_i*(1-e_i);                    % raggio pericentro iniziale
r_p_f = a_f*(1-e_f);                    % raggio pericentro finale

% Range di ricerca per il semiasse maggiore ausiliario
a_aux_min = max(r_p_i, r_p_f) + 100;    
a_aux_max = 200000;                     

% Risoluzione della griglia (Aumenta per più precisione, diminuisci per velocità)
N_a  = 300; 
N_rp = 300; 

a_aux_vec = linspace(a_aux_min, a_aux_max, N_a);

% Variabili per memorizzare l'ottimo
dv_min = Inf;
opt_a_aux = 0;
opt_e_aux = 0;

disp('Scansione 2D (a_aux, r_p_aux) in corso. Attendere...');
tic; % Avvio cronometro

for k = 1:N_a
    a_aux_k = a_aux_vec(k);
    
    % Il raggio di pericentro ausiliario deve essere:
    % 1) Maggiore del raggio terrestre + atmosfera (es. 6500 km)
    % 2) Minore o uguale ad a_aux_k (per avere e_aux tra 0 e 1)
    r_p_aux_vec = linspace(6500, a_aux_k, N_rp);
    
    for j = 1:N_rp
        r_p_aux_k = r_p_aux_vec(j);
        
        % Calcolo eccentricità
        e_aux_k = 1 - r_p_aux_k / a_aux_k;
        
        % Se per qualche errore numerico esce un'iperbole o un cerchio perfetto anomalo, salta
        if e_aux_k >= 1 || e_aux_k < 0
            continue;
        end
        
        % M2: Bitangente 'pa' (da r_p_i ad apocentro ausiliario)
        [dv2_1, dv2_2, ~] = bitangentTransfer(a_i, e_i, a_aux_k, e_aux_k, 'pa', mu);
        
        % M3: Cambio piano
        [dv3, om_tmp, ~] = changeOrbitalPlane(a_aux_k, e_aux_k, i_i, OM_i, om_i, i_f, OM_f, mu);
        
        % Check preliminare M4 per evitare orbite che si incrociano in modo impossibile per 'ap'
        r_a_aux = a_aux_k*(1+e_aux_k);
        if r_a_aux < r_p_f
            continue;
        end
        
        % M4: Bitangente 'ap' (da r_a_aux a r_p_f)
        [dv4_1, dv4_2, ~] = bitangentTransfer(a_aux_k, e_aux_k, a_f, e_f, 'ap', mu);
        
        % M5: Rotazione pericentro
        [dv5, ~, ~] = changePericenterArg(a_f, e_f, om_tmp, om_f, mu);
        
        % Somma totale costi
        dv_tot_current = abs(dv2_1) + abs(dv2_2) + abs(dv3) + abs(dv4_1) + abs(dv4_2) + abs(dv5);
        
        % Check minimo
        if dv_tot_current < dv_min
            dv_min = dv_tot_current;
            opt_a_aux = a_aux_k;
            opt_e_aux = e_aux_k;
        end
    end
end
tempo_calc = toc;

a_aux = opt_a_aux;
e_aux = opt_e_aux;

fprintf('\n=== OTTIMIZZAZIONE COMPLETATA IN %.2f SECONDI ===\n', tempo_calc);
fprintf('  a_aux ottimale = %.4f km\n', a_aux);
fprintf('  e_aux ottimale = %.6f\n', e_aux);
fprintf('  DeltaV TOTALE  = %.6f km/s\n', dv_min);

% -------------------------------------------------------------------------
% 3. RICALCOLO PARAMETRI ESATTI PER IL PLOT
% -------------------------------------------------------------------------
[~, ~, dt_m2] = bitangentTransfer(a_i, e_i, a_aux, e_aux, 'pa', mu);
[~, om_after_plane, th_plane] = changeOrbitalPlane(a_aux, e_aux, i_i, OM_i, om_i, i_f, OM_f, mu);
[~, ~, dt_m4] = bitangentTransfer(a_aux, e_aux, a_f, e_f, 'ap', mu);

% -------------------------------------------------------------------------
% 4. PLOT 3D COMPLETO
% -------------------------------------------------------------------------
figure('Name','Trasferimento Orbitale - Soluzione Ellittica','Color','k','Position',[150 150 1200 900]);
ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal; view(35, 25);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Sequenza Manovre (Orbita Ausiliaria Ellittica)','Color','w','FontSize',14);

% Terra
[XE, YE, ZE] = sphere(50);
surf(XE*6371, YE*6371, ZE*6371, 'EdgeColor', 'none', 'FaceColor', [0.2 0.5 1], 'FaceAlpha', 0.6);

dth = deg2rad(1);

% 1. Orbita Iniziale
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [0.2 0.6 1], '-', 1.5, '1. Orbita Iniziale');

% 2. Trasferimento Salita (Bitangente 'pa')
% Da pericentro iniziale (th=0) ad apocentro ausiliario (th=pi)
r_a_aux_calc = a_aux*(1+e_aux);
a_t1 = (r_p_i + r_a_aux_calc)/2;
e_t1 = (r_a_aux_calc - r_p_i)/(r_a_aux_calc + r_p_i);
plotOrbit(a_t1, e_t1, i_i, OM_i, om_i, 0, pi, dth, mu, [1 0.6 0.1], '-', 2, '2. Trasferimento Salita');

% 3. Orbita Ausiliaria (pre-piano)
plotOrbit(a_aux, e_aux, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [1 1 0], '--', 1.5, '3. Orbita Aux (Pre-Piano)');

% 4. Orbita Ausiliaria (post-piano)
plotOrbit(a_aux, e_aux, i_f, OM_f, om_after_plane, 0, 2*pi, dth, mu, [0.8 0.2 0.8], '-.', 1.5, '4. Orbita Aux (Post-Piano)');

% 5. Trasferimento Discesa (Bitangente 'ap')
% Da apocentro ausiliario (th=pi) a pericentro finale (th=2*pi)
a_t2 = (r_a_aux_calc + r_p_f)/2;
e_t2 = (r_a_aux_calc - r_p_f)/(r_a_aux_calc + r_p_f);
plotOrbit(a_t2, e_t2, i_f, OM_f, om_after_plane, pi, 2*pi, dth, mu, [1 0.3 0.3], ':', 2, '5. Trasferimento Discesa');

% 6. Orbita Finale (Pre-Rotazione Pericentro)
plotOrbit(a_f, e_f, i_f, OM_f, om_after_plane, 0, 2*pi, dth, mu, [0.4 1 0.4], '--', 1, '6. Orbita Finale (Pre-Rotazione)');

% 7. Orbita Finale Target
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, dth, mu, [0 1 0], '-', 2, '7. Orbita Target');

legend('show','Location','eastoutside','TextColor','w','Color','k','FontSize',11);
hold off;

% =========================================================================
% FUNZIONI DI SUPPORTO (assicurati di usare le tue per i calcoli esatti)
% =========================================================================

function plotOrbit(a, e, i, OM, om, th_start, th_end, dth, mu, color, linestyle, linewidth, name)
    th = th_start:dth:th_end;
    if th(end) ~= th_end
        th = [th, th_end];
    end
    r = (a*(1 - e^2)) ./ (1 + e*cos(th));
    r_pqw = [r.*cos(th); r.*sin(th); zeros(1, length(th))];
    R3_OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
    R1_i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
    T_pqw2ijk = R3_OM * R1_i * R3_om;
    r_ijk = T_pqw2ijk * r_pqw;
    plot3(r_ijk(1,:), r_ijk(2,:), r_ijk(3,:), 'Color', color, 'LineStyle', linestyle, ...
          'LineWidth', linewidth, 'DisplayName', name);
end