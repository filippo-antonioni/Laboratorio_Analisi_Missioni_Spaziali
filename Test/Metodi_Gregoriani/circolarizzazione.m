% =========================================================================
%   TRASFERIMENTO ORBITALE - SEQUENZA 5 MANOVRE 
%   OPZIONE: ORBITA AUSILIARIA CIRCOLARE (e_aux = 0)
% =========================================================================
clear; 
clc; 
close all;

% -------------------------------------------------------------------------
% 1. CARICAMENTO DATI
% -------------------------------------------------------------------------
% NOTA: Sostituisci il blocco try-catch con i tuoi dati reali se non hai
% il file a disposizione per il test.
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
% 2. OTTIMIZZAZIONE ORBITA AUSILIARIA CIRCOLARE (e_aux = 0)
% -------------------------------------------------------------------------
r_p_i = a_i*(1-e_i);                    % raggio pericentro orbita iniziale
r_p_f = a_f*(1-e_f);                    % raggio pericentro orbita finale

% L'orbita ausiliaria deve essere più grande dei pericentri per avere senso
% nei trasferimenti bitangenti 'pa' e 'ap'.
a_aux_min = max(r_p_i, r_p_f) + 500;    
a_aux_max = 200000;                     
N_search  = 2000;
a_aux_vec = linspace(a_aux_min, a_aux_max, N_search);
dv_total_vec = zeros(1, N_search);

disp('Scansione di a_aux (circolare) in corso...');

for k = 1:N_search
    a_aux_k = a_aux_vec(k);
    e_aux_k = 0; % FORZIAMO L'ORBITA CIRCOLARE
    
    % Manovra 2: bitangente 'pa' (da pericentro iniziale ad apocentro/orbita circolare)
    [dv2_1, dv2_2, ~] = bitangentTransfer(a_i, e_i, a_aux_k, e_aux_k, 'pa', mu);
    
    % Manovra 3: cambio piano su orbita circolare
    [dv3, om_tmp, th_plane_tmp] = changeOrbitalPlane(a_aux_k, e_aux_k, i_i, OM_i, om_i, i_f, OM_f, mu);
    
    % Manovra 4: bitangente 'ap' (da apocentro/orbita circolare a pericentro finale)
    [dv4_1, dv4_2, ~] = bitangentTransfer(a_aux_k, e_aux_k, a_f, e_f, 'ap', mu);
    
    % Manovra 5: Rotazione pericentro
    [dv5, ~, ~] = changePericenterArg(a_f, e_f, om_tmp, om_f, mu);
    
    % COSTO TOTALE
    dv_total_vec(k) = abs(dv2_1) + abs(dv2_2) + abs(dv3) + abs(dv4_1) + abs(dv4_2) + abs(dv5);
end

% Trova il minimo
[dv_min, idx_opt] = min(dv_total_vec);
a_aux = a_aux_vec(idx_opt);
e_aux = 0; % Circolare

fprintf('\n=== OTTIMIZZAZIONE COMPLETATA ===\n');
fprintf('  a_aux ottimale = %.4f km\n', a_aux);
fprintf('  e_aux          = %.6f (Circolare)\n', e_aux);
fprintf('  DeltaV TOTALE  = %.6f km/s\n', dv_min);

% -------------------------------------------------------------------------
% 3. CALCOLO DETTAGLIATO
% -------------------------------------------------------------------------
% Ricalcoliamo i parametri esatti per l'orbita ottimale
[dv2_1, dv2_2, dt_m2] = bitangentTransfer(a_i, e_i, a_aux, e_aux, 'pa', mu);
[dv3, om_after_plane, th_plane] = changeOrbitalPlane(a_aux, e_aux, i_i, OM_i, om_i, i_f, OM_f, mu);
[dv4_1, dv4_2, dt_m4] = bitangentTransfer(a_aux, e_aux, a_f, e_f, 'ap', mu);
[dv5, thi_vec, thf_vec] = changePericenterArg(a_f, e_f, om_after_plane, om_f, mu);

% -------------------------------------------------------------------------
% 4. PLOT 3D COMPLETO
% -------------------------------------------------------------------------
figure('Name','Trasferimento Orbitale - Soluzione Circolare','Color','k','Position',[100 100 1200 900]);
ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal; view(35, 25);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Sequenza di Manovre Orbitali (Orbita Ausiliaria Circolare)','Color','w','FontSize',14);

% Terra
[XE, YE, ZE] = sphere(50);
surf(XE*6371, YE*6371, ZE*6371, 'EdgeColor', 'none', 'FaceColor', [0.2 0.5 1], 'FaceAlpha', 0.6);

dth = deg2rad(1);

% 1. Orbita Iniziale
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [0.2 0.6 1], '-', 1.5, '1. Orbita Iniziale');

% 2. Trasferimento Salita (Bitangente 'pa')
% Da pericentro iniziale (th=0) ad apocentro di trasferimento (th=pi)
r_a_t1 = a_aux; % Essendo circolare, ra_aux = a_aux
a_t1 = (r_p_i + r_a_t1)/2;
e_t1 = (r_a_t1 - r_p_i)/(r_a_t1 + r_p_i);
plotOrbit(a_t1, e_t1, i_i, OM_i, om_i, 0, pi, dth, mu, [1 0.6 0.1], '-', 2, '2. Trasferimento Salita');

% 3. Orbita Ausiliaria (pre-piano)
plotOrbit(a_aux, e_aux, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [1 1 0], '--', 1.5, '3. Orbita Aux (Pre-Piano)');

% 4. Orbita Ausiliaria (post-piano)
plotOrbit(a_aux, e_aux, i_f, OM_f, om_after_plane, 0, 2*pi, dth, mu, [0.8 0.2 0.8], '-.', 1.5, '4. Orbita Aux (Post-Piano)');

% 5. Trasferimento Discesa (Bitangente 'ap')
% Da apocentro di trasferimento (th=pi) a pericentro finale (th=2*pi o 0)
a_t2 = (a_aux + r_p_f)/2;
e_t2 = (a_aux - r_p_f)/(a_aux + r_p_f);
% L'orbita scende da r_a_t2 a r_p_t2, quindi l'anomalia va da pi a 2*pi
plotOrbit(a_t2, e_t2, i_f, OM_f, om_after_plane, pi, 2*pi, dth, mu, [1 0.3 0.3], ':', 2, '5. Trasferimento Discesa');

% 6. Orbita Finale (Pre-Rotazione Pericentro)
plotOrbit(a_f, e_f, i_f, OM_f, om_after_plane, 0, 2*pi, dth, mu, [0.4 1 0.4], '--', 1, '6. Orbita Finale (Pre-Rotazione)');

% 7. Orbita Finale Target
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, dth, mu, [0 1 0], '-', 2, '7. Orbita Target');

legend('show','Location','eastoutside','TextColor','w','Color','k','FontSize',11);
hold off;


% =========================================================================
% FUNZIONI DI SUPPORTO (da inserire alla fine del file o salvare a parte)
% =========================================================================

function plotOrbit(a, e, i, OM, om, th_start, th_end, dth, mu, color, linestyle, linewidth, name)
    % PLOTORBIT: Traccia un arco di orbita
    th = th_start:dth:th_end;
    if th(end) ~= th_end
        th = [th, th_end];
    end
    
    r = (a*(1 - e^2)) ./ (1 + e*cos(th));
    
    % Perifocali
    r_pqw = [r.*cos(th); r.*sin(th); zeros(1, length(th))];
    
    % Matrici di rotazione
    R3_OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
    R1_i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
    
    T_pqw2ijk = R3_OM * R1_i * R3_om;
    
    r_ijk = T_pqw2ijk * r_pqw;
    
    plot3(r_ijk(1,:), r_ijk(2,:), r_ijk(3,:), 'Color', color, 'LineStyle', linestyle, ...
          'LineWidth', linewidth, 'DisplayName', name);
end