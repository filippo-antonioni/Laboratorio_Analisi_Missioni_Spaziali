% =========================================================================
%   TRASFERIMENTO ORBITALE - SEQUENZA 5 MANOVRE (OPZIONE 1 FIXED)
%   1) Partenza dal punto assegnato dell'orbita iniziale
%   2) Bitangente 'pa' verso orbita ausiliaria ampia (su a_aux)
%   3) Cambio di piano sull'orbita ausiliaria (vicino all'apocentro)
%   4) Bitangente 'ap' per adattare le dimensioni all'orbita finale
%   5) Rotazione dell'argomento del pericentro
% =========================================================================
clear; 
clc; 
close all;

% -------------------------------------------------------------------------
% 1. CARICAMENTO DATI
% -------------------------------------------------------------------------
T  = load("DatiSC1-2026.txt");
mu = 398600;  % [km^3/s^2]
gruppo = 27;

% --- Orbita target / finale ---
rr_f   = T(gruppo, 8:10)';
vv_f   = T(gruppo, 11:13)';
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

% --- Orbita iniziale ---
orbita_iv = T(gruppo,2:7);
a_i  = orbita_iv(1);
e_i  = orbita_iv(2);
i_i  = orbita_iv(3);
OM_i = orbita_iv(4);
om_i = orbita_iv(5);
th_i = orbita_iv(6);

% -------------------------------------------------------------------------
% 2. OTTIMIZZAZIONE ORBITA AUSILIARIA (ciclo for su a_aux)
% -------------------------------------------------------------------------
r_p_i = a_i*(1-e_i);                    % raggio pericentro orbita iniziale
r_a_f = a_f*(1+e_f);                    % raggio apocentro orbita finale

a_aux_min = (r_p_i + r_a_f) / 2;        % la più piccola possibile 
a_aux_max = 800000;                     % a_aux massima (apocentro arriverà a ~1.6M km)
N_search  = 5000;
a_aux_vec = linspace(a_aux_min, a_aux_max, N_search);
dv_total_vec = zeros(1, N_search);

disp('Scansione di a_aux in corso...');
for k = 1:N_search
    a_aux_k = a_aux_vec(k);
    e_aux_k = 1 - r_p_i / a_aux_k;     
    
    % Controllo orbite aperte o impossibili
    if e_aux_k >= 1 || e_aux_k < 0
        dv_total_vec(k) = Inf;
        continue;
    end
    
    % Manovra 2: bitangente 'pa'
    [dv2_1, dv2_2, ~] = bitangentTransfer(a_i, e_i, a_aux_k, e_aux_k, 'pa', mu);
    
    % Manovra 3: cambio piano 
    % Rimosso il try-catch per performance. Le funzioni DEVONO gestire bene gli input
    [dv3, om_tmp, th_plane_tmp] = changeOrbitalPlane(a_aux_k, e_aux_k, i_i, OM_i, om_i, i_f, OM_f, mu);
    
    % Manovra 4: bitangente 'ap'
    r_a_aux = a_aux_k*(1+e_aux_k);
    r_p_f_test = a_f*(1-e_f);
    if r_a_aux < r_p_f_test
        dv_total_vec(k) = Inf;
        continue;
    end
    [dv4_1, dv4_2, ~] = bitangentTransfer(a_aux_k, e_aux_k, a_f, e_f, 'ap', mu);
    
    % Manovra 5: Rotazione pericentro (aggiunta nell'ottimizzazione per equità con Opz 2)
    [dv5, ~, ~] = changePericenterArg(a_f, e_f, om_tmp, om_f, mu);
    
    % COSTO TOTALE (Rigido uso di abs!)
    dv_total_vec(k) = abs(dv2_1) + abs(dv2_2) + abs(dv3) + abs(dv4_1) + abs(dv4_2) + abs(dv5);
end

% Trova il minimo
[dv_min, idx_opt] = min(dv_total_vec);
a_aux = a_aux_vec(idx_opt);
e_aux = 1 - r_p_i / a_aux;

fprintf('\n=== OTTIMIZZAZIONE COMPLETATA ===\n');
fprintf('  a_aux ottimale = %.4f km\n', a_aux);
fprintf('  e_aux          = %.6f\n',    e_aux);
fprintf('  DeltaV TOTALE  = %.6f km/s\n', dv_min);

% -------------------------------------------------------------------------
% 3. CALCOLO DETTAGLIATO E TEMPI DI VOLO (TOF) CORRETTI
% -------------------------------------------------------------------------
% Ricalcoliamo i parametri esatti per l'orbita ottimale
[dv2_1, dv2_2, dt_m2] = bitangentTransfer(a_i, e_i, a_aux, e_aux, 'pa', mu);
[dv3, om_after_plane, th_plane] = changeOrbitalPlane(a_aux, e_aux, i_i, OM_i, om_i, i_f, OM_f, mu);
[dv4_1, dv4_2, dt_m4] = bitangentTransfer(a_aux, e_aux, a_f, e_f, 'ap', mu);
[dv5, thi_vec, thf_vec] = changePericenterArg(a_f, e_f, om_after_plane, om_f, mu);

% Scegliamo la rotazione vicino all'apocentro (minimo dV reale)
[~, idx5] = min(abs(thi_vec - pi));
th_manov5 = thi_vec(idx5);
th_after5 = thf_vec(idx5);

% --- TIMELINE DEI TOF (La vera correzione fisica) ---
% 1. Coasting iniziale (da posizione attuale al pericentro per partire)
dt_coast1 = TOF(a_i, e_i, th_i, 0, mu);

% 2. Trasferimento verso apocentro orbita aux (dt_m2)

% 3. Coasting su orbita aux: dall'apocentro (pi) fino al nodo (th_plane)
dt_attesa_nodo = TOF(a_aux, e_aux, pi, th_plane, mu);

% 4. Coasting su orbita aux: dal nodo (th_plane) di nuovo all'apocentro (pi) per scendere
dt_attesa_discesa = TOF(a_aux, e_aux, th_plane, pi, mu);

% 5. Trasferimento verso orbita finale (dt_m4)

% 6. Coasting su orbita finale: dal pericentro (0) fino al punto di manovra 5 (th_manov5)
dt_attesa_rotazione = TOF(a_f, e_f, 0, th_manov5, mu);

% 7. Coasting finale: dal punto post-manovra 5 fino alla posizione target (th_f)
dt_coast_finale = TOF(a_f, e_f, th_after5, th_f, mu);

% Somma tempi totale
dt_TOT = dt_coast1 + dt_m2 + dt_attesa_nodo + dt_attesa_discesa + dt_m4 + dt_attesa_rotazione + dt_coast_finale;

fprintf('\n========================================================\n');
fprintf('   BUDGET TEMPI E DELTA-V\n');
fprintf('========================================================\n');
fprintf('  Tempo Coasting iniziale:   %.2f giorni\n', dt_coast1/86400);
fprintf('  Tempo Salita (M2):         %.2f giorni\n', dt_m2/86400);
fprintf('  Tempo Attesa Nodo:         %.2f giorni\n', dt_attesa_nodo/86400);
fprintf('  Tempo Attesa Discesa:      %.2f giorni\n', dt_attesa_discesa/86400);
fprintf('  Tempo Discesa (M4):        %.2f giorni\n', dt_m4/86400);
fprintf('  Tempo Attesa Rotazione:    %.2f giorni\n', dt_attesa_rotazione/86400);
fprintf('  Tempo Coasting finale:     %.2f giorni\n', dt_coast_finale/86400);
fprintf('  -----------------------------------------------\n');
fprintf('  TEMPO DI VOLO TOTALE:      %.2f giorni\n', dt_TOT/86400);
fprintf('========================================================\n');

% =========================================================================
% 4. PLOT 3D (Ora con l'orbita di discesa!)
% =========================================================================
figure('Name','Trasferimento Orbitale - Vista 3D (Fixed)','Color','k','Position',[100 100 1200 900]);
ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal; view(35, 25);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Sequenza di Manovre Orbitali','Color','w','FontSize',14);
dth = deg2rad(0.5);

% Orbita Iniziale
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi-dth, dth, mu);
h_init = findobj(gca,'Type','Line'); set(h_init(1),'Color',[0.2 0.6 1],'DisplayName','Orbita Iniziale');

% Trasferimento Salita
r_a_aux2 = a_aux*(1+e_aux);
a_t2 = (r_p_i + r_a_aux2)/2;
e_t2 = (r_a_aux2 - r_p_i)/(r_a_aux2 + r_p_i);
plotOrbit(a_t2, e_t2, i_i, OM_i, om_i, 0, pi, dth, mu);
h_t2 = findobj(gca,'Type','Line'); set(h_t2(1),'Color',[1 0.6 0.1],'DisplayName','Trasferimento Salita');

% Orbita Ausiliaria (pre-piano)
plotOrbit(a_aux, e_aux, i_i, OM_i, om_i, 0, 2*pi-dth, dth, mu);
h_aux1 = findobj(gca,'Type','Line'); set(h_aux1(1),'Color',[1 1 0],'LineStyle','--','DisplayName','Orbita Aux (pre-piano)');

% Orbita Ausiliaria (post-piano)
plotOrbit(a_aux, e_aux, i_f, OM_f, om_after_plane, 0, 2*pi-dth, dth, mu);
h_aux2 = findobj(gca,'Type','Line'); set(h_aux2(1),'Color',[0.8 0.2 0.8],'LineStyle','-.','DisplayName','Orbita Aux (post-piano)');

% Trasferimento Discesa (Aggiunto!)
r_p_f2 = a_f*(1-e_f);
a_t4 = (r_a_aux2 + r_p_f2)/2;
e_t4 = (r_a_aux2 - r_p_f2)/(r_a_aux2 + r_p_f2);
plotOrbit(a_t4, e_t4, i_f, OM_f, om_after_plane, pi, 2*pi, dth, mu);
h_t4 = findobj(gca,'Type','Line'); set(h_t4(1),'Color',[1 0.3 0.3],'LineStyle',':','LineWidth',2,'DisplayName','Trasferimento Discesa');

% Orbita Finale
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi-dth, dth, mu);
h_fin = findobj(gca,'Type','Line'); set(h_fin(1),'Color',[0.1 1 0.1],'LineWidth', 2, 'DisplayName','Orbita Finale Target');

legend('show','Location','bestoutside','TextColor','w','Color','k','FontSize',10);