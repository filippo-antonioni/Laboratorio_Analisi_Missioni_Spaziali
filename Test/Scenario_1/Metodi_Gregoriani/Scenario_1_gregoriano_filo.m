%% Opzione 2: Ottimizzazione Trasferimento con Orbita di Appoggio Estesa
clear
close all
clc
gruppo=23;
disp('--- AVVIO OPZIONE 2: Ricerca orbita di trasferimento AMPIA (Grid Search) ---');

% --- 1. CARICAMENTO DATI ---
T = load("DatiSC1-2026.txt");
mu = 398600;

% Orbita target/finale
rr_f = T(gruppo,8:10)'; 
vv_f = T(gruppo,11:13)';
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

% Orbita iniziale
orbita_iv = T(gruppo,2:7);
a_i = orbita_iv(1); 
e_i = orbita_iv(2); 
i_i = orbita_iv(3);
OM_i = orbita_iv(4); 
om_i = orbita_iv(5);
th_i = orbita_iv(6);

% --- 2. RICERCA DEL MINIMO GLOBALE (Ciclo For) ---
r_p_i = a_i * (1 - e_i); % Raggio pericentro iniziale
r_a_i = a_i * (1 + e_i); % Raggio apocentro iniziale
r_p_f = a_f * (1 - e_f); % Raggio pericentro finale

% Nell'Opzione 1 avevi a_max = 800000. L'apocentro corrispondente è 2*800000 - r_p_i
r_a_max = 2 * 800000 - r_p_i; 
r_a_test_vec = linspace(max(r_a_i, r_p_f), r_a_max, 1500000);

DeltaV_minimo = inf;
fprintf('Scansione di %d orbite di appoggio in corso...\n', length(r_a_test_vec));

for idx = 1:length(r_a_test_vec)
    r_a_test = r_a_test_vec(idx);
    
    % Definiamo l'orbita intermedia (tangente al pericentro iniziale)
    a_int = (r_p_i + r_a_test) / 2;
    e_int = (r_a_test - r_p_i) / (r_a_test + r_p_i);
    
    if e_int >= 1 || e_int < 0
        continue; % Scarta le orbite non chiuse
    end
    
    
    % FASE 1: Bitangente 1 ('pa' - da Pericentro i ad Apocentro int)
    [dV1_A, dV1_B, dt1] = bitangentTransfer(a_i, e_i, a_int, e_int, 'pa', mu);
    costo_bitang1 = abs(dV1_A) + abs(dV1_B);
    
    % FASE 2: Cambio Piano (sull'orbita intermedia)
    [dV_plane, om_fun, theta_plane] = changeOrbitalPlane(a_int, e_int, i_i, OM_i, om_i, i_f, OM_f, mu);
    
    % FASE 3: Bitangente 2 ('ap' - da Apocentro int a Pericentro f)
    [dV2_A, dV2_B, dt2] = bitangentTransfer(a_int, e_int, a_f, e_f, 'ap', mu);
    costo_bitang2 = abs(dV2_A) + abs(dV2_B);
    
    % FASE 4: Cambio Anomalia Pericentro (sull'orbita finale)
    [dV_arg, thi_fun, thf_fun] = changePericenterArg(a_f, e_f, om_fun, om_f, mu);
    
    % COSTO TOTALE
    costo_totale = costo_bitang1 + abs(dV_plane) + costo_bitang2 + abs(dV_arg); % ho aggiunto gli abs - simone
    
    % Aggiornamento del minimo
    if costo_totale < DeltaV_minimo
        DeltaV_minimo = costo_totale;
        
        best_a_int = a_int; 
        best_e_int = e_int;
        best_dv1 = costo_bitang1;
        best_dv_plane = dV_plane;
        best_dv2 = costo_bitang2;
        best_dv_arg = dV_arg;
        
        best_dt1 = dt1;
        best_dt2 = dt2;
        best_om_fun = om_fun; 
        best_theta_plane = theta_plane;
        [best_thi_fun,index_best_thi_fun] = min(thi_fun);
        best_thf_fun = thf_fun(index_best_thi_fun);
    end
end

% --- 3. CALCOLO DEI TEMPI DI VOLO (TOF) ---
deltat_coast1 = TOF(a_i, e_i, th_i, 0, mu); % Coasting fino al pericentro iniziale
deltat_attesa1 = TOF(best_a_int, best_e_int, pi, best_theta_plane, mu); % Da apocentro a nodo
deltat_attesa2 = TOF(best_a_int, best_e_int, best_theta_plane, pi, mu); % Da nodo di nuovo ad apocentro
deltat_attesa3 = TOF(a_f, e_f, 0, best_thi_fun, mu); % Da arrivo al pericentro fino al cambio om
deltat_coast_finale = TOF(a_f, e_f, best_thf_fun, th_f, mu); % Fino alla posizione finale

TOF_Totale = deltat_coast1 + best_dt1 + deltat_attesa1 + deltat_attesa2 + best_dt2 + deltat_attesa3 + deltat_coast_finale;
TOF_giorni = TOF_Totale / (24*3600);

% --- STAMPE A SCHERMO ---
fprintf('\n--- RISULTATI OPZIONE 2 (Ottimo con Orbita di Appoggio) ---\n');
fprintf('Semiasse Orbita Appoggio: %.2f km\n', best_a_int);
fprintf('Eccentricità Appoggio:    %.4f\n\n', best_e_int);

fprintf('DeltaV Bitangente 1:      %.4f km/s\n', best_dv1);
fprintf('DeltaV Cambio Piano:      %.4f km/s\n', best_dv_plane);
fprintf('DeltaV Bitangente 2:      %.4f km/s\n', best_dv2);
fprintf('DeltaV Cambio Pericen.:   %.4f km/s\n', best_dv_arg);
fprintf('DELTAV TOTALE:            %.4f km/s\n\n', DeltaV_minimo);

fprintf('--- TEMPI DI VOLO ---\n');
fprintf('Coasting Iniziale:        %.2f giorni\n', deltat_coast1/(24*3600));
fprintf('Trasferimento 1 (Salita): %.2f giorni\n', best_dt1/(24*3600));
fprintf('Attesa (verso il nodo):   %.2f giorni\n', deltat_attesa1/(24*3600));
fprintf('Attesa (verso discesa):   %.2f giorni\n', deltat_attesa2/(24*3600));
fprintf('Trasferimento 2 (Disce.): %.2f giorni\n', best_dt2/(24*3600));
fprintf('Attesa (verso manovra):   %.2f giorni\n', deltat_attesa3/(24*3600));
fprintf('Coasting Finale:          %.2f giorni\n', deltat_coast_finale/(24*3600));
fprintf('TEMPO DI VOLO TOTALE:     %.2f giorni\n\n', TOF_giorni);

% =========================================================================
% --- 4. GRAFICA AVANZATA 3D (Tutte le Orbite e le Manovre) ---
% =========================================================================
figure('Name', 'Opzione 2: Sequenza Completa con Tutte le Orbite', 'NumberTitle', 'off');
hold on; grid on; axis equal; view(3);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

th_vec = linspace(0, 2*pi, 300);
get_c = @(a,e,i,OM,om,th,idx) subsref(par2car(a,e,i,OM,om,th,mu), struct('type','()','subs',{{idx}}));

% --- CALCOLO DELLA SECONDA ORBITA DI TRASFERIMENTO (Quella mancante!) ---
% Il satellite scende dall'apocentro dell'orbita intermedia al pericentro di quella finale
r_a_int = best_a_int * (1 + best_e_int);
r_p_f = a_f * (1 - e_f);
a_t2 = (r_a_int + r_p_f) / 2;
e_t2 = abs(r_a_int - r_p_f) / (r_a_int + r_p_f);

% --- PLOT DELLE 6 ORBITE ---

% 1. Orbita Iniziale (Blu)
plot3(arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,1), th_vec), ...
      arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,2), th_vec), ...
      arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,3), th_vec), ...
      'b', 'LineWidth', 1.5, 'DisplayName', '1. Orbita Iniziale');

% 2. Appoggio 1: Salita (Arancione tratteggiata - Pre-Cambio Piano)
plot3(arrayfun(@(th) get_c(best_a_int,best_e_int,i_i,OM_i,om_i,th,1), th_vec), ...
      arrayfun(@(th) get_c(best_a_int,best_e_int,i_i,OM_i,om_i,th,2), th_vec), ...
      arrayfun(@(th) get_c(best_a_int,best_e_int,i_i,OM_i,om_i,th,3), th_vec), ...
      'color', [1 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', '2. Appoggio 1 (Piano i)');

% 3. Appoggio 2: Apocentro (Ciano tratto-punto - Post-Cambio Piano)
plot3(arrayfun(@(th) get_c(best_a_int,best_e_int,i_f,OM_f,best_om_fun,th,1), th_vec), ...
      arrayfun(@(th) get_c(best_a_int,best_e_int,i_f,OM_f,best_om_fun,th,2), th_vec), ...
      arrayfun(@(th) get_c(best_a_int,best_e_int,i_f,OM_f,best_om_fun,th,3), th_vec), ...
      'c', 'LineStyle', '-.', 'LineWidth', 1.5, 'DisplayName', '3. Appoggio 2 (Piano f)');

% 4. TRASFERIMENTO DI DISCESA (Rosso puntinato - L'orbita che mancava!)
plot3(arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,best_om_fun,th,1), th_vec), ...
      arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,best_om_fun,th,2), th_vec), ...
      arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,best_om_fun,th,3), th_vec), ...
      'r', 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', '4. Ellisse Discesa (Bitang. 2)');

% 5. Orbita Target Pre-Pericentro (Magenta continua)
plot3(arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,best_om_fun,th,1), th_vec), ...
      arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,best_om_fun,th,2), th_vec), ...
      arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,best_om_fun,th,3), th_vec), ...
      'm', 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', '5. Orbita Target (om sfasato)');

% 6. Orbita Target Finale (Verde spessa)
plot3(arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,1), th_vec), ...
      arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,2), th_vec), ...
      arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,3), th_vec), ...
      'g', 'LineWidth', 2.5, 'DisplayName', '6. Orbita Target Finale');


% --- MARKERS DELLE MANOVRE ---
% M1: Inizio Salita (Pericentro Iniziale)
[r_m1, ~] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu); 
plot3(r_m1(1), r_m1(2), r_m1(3), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'M1: Inizio Salita');

% M2: Cambio Piano (Sull'orbita di appoggio 1)
[r_m2, ~] = par2car(best_a_int, best_e_int, i_i, OM_i, om_i, best_theta_plane, mu); 
plot3(r_m2(1), r_m2(2), r_m2(3), '^', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'DisplayName', 'M2: Cambio Piano');

% M3: Inizio Discesa (Apocentro dell'orbita di appoggio 2)
[r_m3, ~] = par2car(best_a_int, best_e_int, i_f, OM_f, best_om_fun, pi, mu); 
plot3(r_m3(1), r_m3(2), r_m3(3), 's', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'M3: Inizio Discesa');

% M4: Fine Discesa / Circolarizzazione (Pericentro dell'ellisse di discesa)
[r_m4, ~] = par2car(a_t2, e_t2, i_f, OM_f, best_om_fun, 0, mu); 
plot3(r_m4(1), r_m4(2), r_m4(3), 'd', 'MarkerSize', 7, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'DisplayName', 'M4: Fine Discesa');

% M5: Cambio Pericentro (Intersezione tra l'orbita magenta e verde)
[r_m5, ~] = par2car(a_f, e_f, i_f, OM_f, best_om_fun, best_thi_fun, mu); 
plot3(r_m5(1), r_m5(2), r_m5(3), 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'M5: Cambio Pericentro');


% --- TERRA 3D A PROVA DI CRASH ---
R_earth = 6371; 
[xE, yE, zE] = sphere(50); 
try
    load topo topo topomap1; 
    surf(xE * R_earth, yE * R_earth, zE * R_earth, 'FaceColor', 'texturemap', 'CData', topo, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    colormap(topomap1); 
catch
    % Fallback pulito: se topo non carica, crea un globo blu compatto
    surf(xE * R_earth, yE * R_earth, zE * R_earth, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

legend('show', 'Location', 'bestoutside'); 
hold off;