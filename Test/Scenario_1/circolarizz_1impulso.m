clear; 
clc; 
close all;

% --- 1. CARICAMENTO DATI ---
T = load("DatiSC1-2026.txt");
gruppo = 23;
mu = 398600;
R_earth = 6371;

% orbita finale
rr_f = T(gruppo,8:10)'; 
vv_f = T(gruppo,11:13)';
orbita_iv = T(gruppo,2:7);
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

% orbita iniziale
a_i = orbita_iv(1); 
e_i = orbita_iv(2); 
i_i = orbita_iv(3);
OM_i = orbita_iv(4); 
om_i = orbita_iv(5); 
th_i = orbita_iv(6);

% --- 2. SETUP E SEQUENZA MANOVRE ---
r_p_i = a_i * (1 - e_i); 
r_a_i = a_i * (1 + e_i); 
r_p_f = a_f * (1 - e_f); 

% Coasting dal punto di partenza fino all'apocentro
dt_coast1 = TOF(a_i, e_i, th_i, pi, mu);

% Calcolo posizioni e velocità all'apocentro sull'orbita iniziale 
[rr_i_a, vv_i_a] = par2car(a_i, e_i, i_i, OM_i, om_i, pi, mu);

% Calcolo orbita di parcheggio circolare (a_park = r_a_i corretto)
a_park = r_a_i;
e_park = 0;
[rr_park, vv_park] = par2car(a_park, e_park, i_i, OM_i, om_i, pi, mu);
costo_circolarizzazione = norm(vv_i_a - vv_park);

% Calcolo deltav del cambio pian
[dV_plane, om_plane, theta_plane] = changeOrbitalPlane(a_park, e_park, i_i, OM_i, om_i, i_f, OM_f, mu);
costo_piano = abs(dV_plane);
        
% --- RICERCA DEL NODO PIU' VICINO ---
% Trovo i due nodi di intersezione tra i piani
th_nodo1 = mod(theta_plane, 2*pi);
th_nodo2 = mod(theta_plane + pi, 2*pi);

% Scelgo il nodo che incontro prima partendo dalla circolarizzazione (th = pi)
dt_nodo1 = TOF(a_park, e_park, pi, th_nodo1, mu);
dt_nodo2 = TOF(a_park, e_park, pi, th_nodo2, mu);

if dt_nodo1 < dt_nodo2
    th_plane_tmp = th_nodo1;
    dt_coast_park = dt_nodo1;
else
    th_plane_tmp = th_nodo2;
    dt_coast_park = dt_nodo2;
end
        
% --- CAMBIO PERICENTRO SU CIRCOLARE (Istantaneo) ---
[dV_arg, thi_fun, thf_fun] = changePericenterArg(a_park, e_park, om_plane, om_f, mu);

% Essendo un'orbita circolare (e=0), dV_arg = 0. Lo applichiamo istantaneamente 
% nello stesso momento e nello stesso punto spaziale del cambio piano.
dt_coast_periarg = 0; 

% Calcoliamo l'anomalia in cui si trova il satellite nell'istante del cambio 
% piano, ma letta nel NUOVO sistema di riferimento ruotato (con om_f)
th_post_plane_new = mod(th_plane_tmp + om_plane - om_f, 2*pi);
       
% =========================================================================
% OPZIONE 1: DISCESA CON BITANGENTE
% =========================================================================
dv_bitang = [];
[DeltaV1_ap, DeltaV2_ap, Deltat_ap] = bitangentTransfer(a_park, e_park, a_f, e_f, 'ap', mu);
dv_bitang = [dv_bitang; abs(DeltaV1_ap)+abs(DeltaV2_ap)];

[DeltaV1_pa, DeltaV2_pa, Deltat_pa] = bitangentTransfer(a_park, e_park, a_f, e_f, 'pa', mu);
dv_bitang = [dv_bitang; abs(DeltaV1_pa)+abs(DeltaV2_pa)];

[DeltaV1_aa, DeltaV2_aa, Deltat_aa] = bitangentTransfer(a_park, e_park, a_f, e_f, 'aa', mu);
dv_bitang = [dv_bitang; abs(DeltaV1_aa)+abs(DeltaV2_aa)];

[DeltaV1_pp, DeltaV2_pp, Deltat_pp] = bitangentTransfer(a_park, e_park, a_f, e_f, 'pp', mu);
dv_bitang = [dv_bitang; abs(DeltaV1_pp)+abs(DeltaV2_pp)];

[dv_bitang2, idx_bitang2] = min(dv_bitang);

switch idx_bitang2
    case 1
        th_bitang = pi;
        th_arrivo_bitang = 0;
        dt_bitang2 = Deltat_ap;
    case 2
        th_bitang = 0;
        th_arrivo_bitang = pi;
        dt_bitang2 = Deltat_pa;
    case 3
        th_bitang = pi;
        th_arrivo_bitang = pi;
        dt_bitang2 = Deltat_aa;
    case 4
        th_bitang = 0;
        th_arrivo_bitang = 0;
        dt_bitang2 = Deltat_pp;
end

% Coasting dal punto di cambio piano/peri (letto nel nuovo riferimento) fino alla bitangente
dt_coast_park2 = TOF(a_park, e_park, th_post_plane_new, th_bitang, mu);
dt_coast_finale_A = TOF(a_f, e_f, th_arrivo_bitang, th_f, mu);  
  
tempo_tot_A = dt_coast1 + dt_coast_park + dt_coast_periarg + dt_bitang2 + dt_coast_park2 + dt_coast_finale_A;
dV_tot_A = costo_circolarizzazione + costo_piano + abs(dV_arg) + dv_bitang2;

% =========================================================================
% OPZIONE 2: INTERSEZIONE DIRETTA (SE ESISTE)
% =========================================================================
[interseca, th1_int, th2_int] = orbit_intersection(a_park, e_park, om_f, a_f, e_f, om_f);
if interseca 
    if th_post_plane_new < th1_int(1) || th_post_plane_new > th1_int(2)
        dt_coast_intersez = TOF(a_park, e_park, th_post_plane_new, th1_int(1), mu);
        th1_intersez = th1_int(1);
        th2_intersez = th2_int(1);
    else 
        dt_coast_intersez = TOF(a_park, e_park, th_post_plane_new, th1_int(2), mu);
        th1_intersez = th1_int(2);
        th2_intersez = th2_int(2);
    end
        
    [rr_intersez_1, vv_intersez_1] = par2car(a_park, e_park, i_f, OM_f, om_f, th1_intersez, mu);
    [rr_intersez_2, vv_intersez_2] = par2car(a_f, e_f, i_f, OM_f, om_f, th2_intersez, mu); 
    dv_discesa = norm(vv_intersez_1 - vv_intersez_2);
    dt_coast_finale_B = TOF(a_f, e_f, th2_intersez, th_f, mu);
    tempo_tot_B = dt_coast1 + dt_coast_park + dt_coast_periarg + dt_coast_intersez + dt_coast_finale_B;
    dV_tot_B = costo_circolarizzazione + costo_piano + abs(dV_arg) + dv_discesa;
else 
    dV_tot_B = inf; 
    tempo_tot_B = inf;
end 

% =========================================================================
% --- 3. STAMPA RISULTATI DINAMICA E CONFRONTO ---
% =========================================================================
% Determino in automatico la migliore per i plot e il riepilogo
if dV_tot_A <= dV_tot_B
    min_costo = dV_tot_A;
    tof_tot = tempo_tot_A;
    opzione = 'A (Bitangente Ottimizzata)';
else
    min_costo = dV_tot_B;
    tof_tot = tempo_tot_B;
    opzione = 'B (Intersezione Diretta)';
end

fprintf('\n=======================================================\n');
fprintf('   CONFRONTO STRATEGIE SULLA CIRCOLARE\n');
fprintf('=======================================================\n');

% --- STAMPA OPZIONE A ---
fprintf('\n---> OPZIONE A: Discesa con Bitangente\n');
fprintf('DeltaV Totale A: %.4f km/s | TOF: %.2f giorni\n', dV_tot_A, tempo_tot_A/86400);
fprintf('- Impulso 1 (Circolarizzazione):%.4f km/s\n', costo_circolarizzazione);
fprintf('- Impulso 2 (Cambio Piano):     %.4f km/s\n', costo_piano);
fprintf('- Impulso 3 (Cambio Pericentro):%.4f km/s (Matematico)\n', abs(dV_arg));
fprintf('- Impulsi 4+5 (Bitangente):     %.4f km/s\n', dv_bitang2);
fprintf('  [Tempi di Volo A]\n');
fprintf('   Coasting Iniziale:      %.2f gg\n', dt_coast1 / 86400);
fprintf('   Coasting Cambio Piano:  %.2f gg (Verso il nodo più vicino)\n', dt_coast_park / 86400);
fprintf('   Coasting Cambio Peri:   %.2f gg (Istantaneo)\n', dt_coast_periarg / 86400);
fprintf('   Coasting a Bitangente:  %.2f gg\n', dt_coast_park2 / 86400);
fprintf('   Trasferimento Bitang:   %.2f gg\n', dt_bitang2 / 86400);
fprintf('   Coasting Finale:        %.2f gg\n', dt_coast_finale_A / 86400);

% --- STAMPA OPZIONE B ---
if isinf(dV_tot_B)
    fprintf('\n---> OPZIONE B: Intersezione Diretta\n');
    fprintf('NON FATTIBILE: L''orbita di parcheggio non interseca quella target.\n');
else
    dt_coast_intersez_print = tempo_tot_B - dt_coast1 - dt_coast_park - dt_coast_periarg - dt_coast_finale_B;
    fprintf('\n---> OPZIONE B: Intersezione Diretta\n');
    fprintf('DeltaV Totale B: %.4f km/s | TOF: %.2f giorni\n', dV_tot_B, tempo_tot_B/86400);
    fprintf('- Impulso 1 (Circolarizzazione):%.4f km/s\n', costo_circolarizzazione);
    fprintf('- Impulso 2 (Cambio Piano):     %.4f km/s\n', costo_piano);
    fprintf('- Impulso 3 (Cambio Pericentro):%.4f km/s (Matematico)\n', abs(dV_arg));
    fprintf('- Impulso 4 (Salto Diretto):    %.4f km/s\n', dv_discesa);
    fprintf('  [Tempi di Volo B]\n');
    fprintf('   Coasting Iniziale:      %.2f gg\n', dt_coast1 / 86400);
    fprintf('   Coasting Cambio Piano:  %.2f gg\n', dt_coast_park / 86400);
    fprintf('   Coasting Cambio Peri:   %.2f gg (Istantaneo)\n', dt_coast_periarg / 86400);
    fprintf('   Coasting a Intersez:    %.2f gg\n', dt_coast_intersez_print / 86400);
    fprintf('   Coasting Finale:        %.2f gg\n', dt_coast_finale_B / 86400);
end

fprintf('\n=======================================================\n');
fprintf('>>> MIGLIOR STRATEGIA VINCITRICE: %s <<<\n', opzione);
fprintf('=======================================================\n');

% =========================================================================
% --- 4. GRAFICA AVANZATA 3D (Si adatta alla strategia vincente) ---
% =========================================================================
fig = figure('Name', 'Analisi Missione: Trasferimento Misto', 'NumberTitle', 'off', 'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
hold on; grid on; axis equal; view(3); rotate3d on;
xlabel('X [km]', 'FontWeight', 'bold'); ylabel('Y [km]', 'FontWeight', 'bold'); zlabel('Z [km]', 'FontWeight', 'bold');
title(['Traiettoria 3D Ottimizzata - Ha Vinto: ', opzione], 'FontSize', 14, 'FontWeight', 'bold');

th_vec_full = linspace(0, 2*pi, 300);
get_c = @(a,e,i,OM,om,th,idx) subsref(par2car(a,e,i,OM,om,th,mu), struct('type','()','subs',{{idx}}));

% 1. Orbita Iniziale
plot3(arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,3), th_vec_full), ...
      'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', 'Orbita Iniziale');

[r_start,~] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu);
plot3(r_start(1), r_start(2), r_start(3), 'ok', 'MarkerSize', 7, 'MarkerFaceColor', [0 0.4470 0.7410], 'DisplayName', 'Pericentro Iniziale');

% Circolarizzazione all'Apocentro
[r_circ,~] = par2car(a_i, e_i, i_i, OM_i, om_i, pi, mu);
plot3(r_circ(1), r_circ(2), r_circ(3), 'sk', 'MarkerSize', 7, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'DisplayName', 'Circolarizzazione');

% 2. Orbita Appoggio Pre-Piano
plot3(arrayfun(@(th) get_c(a_park,e_park,i_i,OM_i,om_i,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_i,OM_i,om_i,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_i,OM_i,om_i,th,3), th_vec_full), ...
      'Color', [0.9290 0.6940 0.1250], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Pre-Plane)');

% Marker Cambio Piano (Triangolo Grande Sotto)
[r_plane,~] = par2car(a_park, e_park, i_i, OM_i, om_i, th_plane_tmp, mu);
plot3(r_plane(1), r_plane(2), r_plane(3), '^k', 'MarkerSize', 14, 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Cambio Piano');

% 3. Orbita Appoggio Post-Piano
plot3(arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_plane,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_plane,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_plane,th,3), th_vec_full), ...
      'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Post-Plane)');

% Cambio Pericentro (Stella Piccola Sopra, stesso identico punto)
[r_cambio_peri,~] = par2car(a_park, e_park, i_f, OM_f, om_f, th_post_plane_new, mu);
plot3(r_cambio_peri(1), r_cambio_peri(2), r_cambio_peri(3), 'pk', 'MarkerSize', 7, 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'DisplayName', 'Cambio Pericentro (Matematico)');

% 4. Orbita Appoggio Post-Peri
plot3(arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_f,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_f,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_f,th,3), th_vec_full), ...
      'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Post-Peri)');

% --- DISEGNO DINAMICO DELLA MANOVRA FINALE E MARKER DI ARRIVO ---
if opzione(1) == 'A'
    % Ricostruisco l'ellisse della bitangente
    if idx_bitang2 == 1 || idx_bitang2 == 4
        rt2 = a_f*(1-e_f);
    else
        rt2 = a_f*(1+e_f);
    end
    a_t2 = (a_park + rt2)/2;
    e_t2 = abs(a_park - rt2)/(a_park + rt2);
    
    th_vec_t2 = linspace(th_bitang, th_bitang+pi, 150);
    plot3(arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,om_f,th,1), th_vec_t2), ...
          arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,om_f,th,2), th_vec_t2), ...
          arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,om_f,th,3), th_vec_t2), ...
          'Color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Bitangente (Opzione A)');
          
    [r_start_discesa,~] = par2car(a_park, e_park, i_f, OM_f, om_f, th_bitang, mu);
    plot3(r_start_discesa(1), r_start_discesa(2), r_start_discesa(3), 'dk', 'MarkerSize', 7, 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'DisplayName', 'Inizio Discesa (A)');

    % Marker arrivo su Orbita Finale Target
    [r_arrivo,~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_arrivo_bitang, mu);
    plot3(r_arrivo(1), r_arrivo(2), r_arrivo(3), 'p', 'MarkerSize', 12, 'MarkerFaceColor', '#77AC30', 'MarkerEdgeColor', 'w', 'DisplayName', 'Arrivo su Target');
else
    % Se ha vinto l'intersezione diretta
    if exist('th1_intersez', 'var')
        [r_intersez,~] = par2car(a_park, e_park, i_f, OM_f, om_f, th1_intersez, mu);
        plot3(r_intersez(1), r_intersez(2), r_intersez(3), 'hk', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Punto di Intersezione (B)');
        
        % Marker arrivo su Orbita Finale Target (stesso punto)
        [r_arrivo,~] = par2car(a_f, e_f, i_f, OM_f, om_f, th2_intersez, mu);
        plot3(r_arrivo(1), r_arrivo(2), r_arrivo(3), 'p', 'MarkerSize', 12, 'MarkerFaceColor', '#77AC30', 'MarkerEdgeColor', 'w', 'DisplayName', 'Arrivo su Target');
    end
end

% 5. Orbita Finale Target
plot3(arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,3), th_vec_full), ...
      'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.5, 'DisplayName', 'Orbita Finale Target');

% --- TERRA FOTOREALISTICA ---
[xE, yE, zE] = sphere(50); 
try
    load topo topo topomap1;
    surf(xE * R_earth, yE * R_earth, zE * R_earth, 'FaceColor', 'texturemap', 'CData', topo, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    colormap(topomap1);
catch
    surf(xE * R_earth, yE * R_earth, zE * R_earth, 'FaceColor', [0.1 0.4 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
legend('show', 'Location', 'bestoutside', 'FontSize', 10); 
hold off;

%% ANIMAZIONE DELLA MANOVRA
% =========================================================================
% ANIMAZIONE 3D AVANZATA - TRACCIA MULTI-COLORE E ORBITE DI RIFERIMENTO
% =========================================================================
if ~exist('opzione', 'var')
    error('Esegui prima il Codice 2 per caricare i dati nel Workspace!');
end
% --- 1. SETUP SCENA ---
fig_anim = figure('Name', 'Simulazione Dinamica Trasferimento Orbitale', 'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
hold on; grid on; axis equal; view(45, 30);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title({'Missione Gruppo 23: Animazione Sequenza Manovre', ['Strategia Vincente: ', opzione]}, 'FontSize', 14);
% Terra 3D
[xE, yE, zE] = sphere(50);
try
    load topo topo topomap1;
    surf(xE * 6371, yE * 6371, zE * 6371, 'FaceColor', 'texturemap', 'CData', topo, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    colormap(topomap1);
catch
    surf(xE * 6371, yE * 6371, zE * 6371, 'FaceColor', [0.1 0.4 0.8], 'EdgeColor', 'none');
end
% --- 2. PLOT ORBITE DI RIFERIMENTO (Tratteggiate) ---
dth_plot = linspace(0, 2*pi, 300);
% Orbita Iniziale
plotOrbitStatic(a_i, e_i, i_i, OM_i, om_i, mu, [0 0.45 0.74], '--', 1, 'Rif. Iniziale');
% Orbita Parcheggio Pre-Piano
plotOrbitStatic(a_park, e_park, i_i, OM_i, om_i, mu, [0.93 0.69 0.13], '--', 1, 'Rif. Parcheggio Pre-Piano');
% Orbita Parcheggio Post-Piano
plotOrbitStatic(a_park, e_park, i_f, OM_f, om_plane, mu, [0.49 0.18 0.56], '--', 1, 'Rif. Parcheggio Post-Piano');
% Orbita Finale Target
plotOrbitStatic(a_f, e_f, i_f, OM_f, om_f, mu, [0.47 0.67 0.19], '--', 1, 'Rif. Finale');
% --- 3. DEFINIZIONE SEGMENTI E COLORI TRACCIA ---
% Ogni riga: {a, e, i, OM, om, th_start, th_end, nome, colore_traccia}
segmenti = {};
% Seg 1: Da th_i all'apocentro (Orbita Iniziale)
th_e = pi; if th_e < th_i; th_e = th_e + 2*pi; end
segmenti{1} = {a_i, e_i, i_i, OM_i, om_i, th_i, th_e, 'Coasting Iniziale', [0 0.44 0.74]};
% Seg 2: Dal pi al nodo di cambio piano (Orbita Parcheggio 1)
th_s = pi; th_e = th_plane_tmp; if th_e < th_s; th_e = th_e + 2*pi; end
segmenti{2} = {a_park, e_park, i_i, OM_i, om_i, th_s, th_e, 'Verso Nodo Cambio Piano', [0.85 0.32 0.1]};
% Seg 3: Dal nodo al punto di inizio discesa (Orbita Parcheggio 2 - Post Piano)
th_s = th_post_plane_new; 
if opzione(1) == 'A'
    th_e = th_bitang;
else
    th_e = th1_intersez;
end
if th_e < th_s; th_e = th_e + 2*pi; end
segmenti{3} = {a_park, e_park, i_f, OM_f, om_f, th_s, th_e, 'Allineamento Discesa', [0.49 0.18 0.56]};
% Seg 4: La Manovra Finale (Bitangente o Intersezione)
if opzione(1) == 'A'
    segmenti{4} = {a_t2, e_t2, i_f, OM_f, om_f, th_bitang, th_bitang+pi, 'Trasferimento Bitangente', [0.63 0.07 0.18]};
    th_fine_manovra = th_arrivo_bitang;
else
    segmenti{4} = {a_park, e_park, i_f, OM_f, om_f, th1_intersez, th1_intersez+0.01, 'Salto Intersezione', [0 0 0]};
    th_fine_manovra = th2_intersez;
end
% Seg 5: Coasting finale fino a th_f (Orbita Target)
th_s = th_fine_manovra; th_e = th_f; if th_e < th_s; th_e = th_e + 2*pi; end
segmenti{5} = {a_f, e_f, i_f, OM_f, om_f, th_s, th_e, 'Arrivo a Target', [0.46 0.67 0.18]};
% --- 4. CICLO DI ANIMAZIONE ---
h_sat = plot3(NaN, NaN, NaN, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Satellite');
step = 100; % Punti per ogni segmento
for s = 1:length(segmenti)
    seg = segmenti{s};
    th_v = linspace(seg{6}, seg{7}, step);
    
    % Crea una nuova linea per la traccia di questo segmento (cambio colore)
    h_trail = plot3(NaN, NaN, NaN, 'Color', seg{9}, 'LineWidth', 2, 'DisplayName', seg{8});
    t_x = []; t_y = []; t_z = [];
    
    for k = 1:length(th_v)
        r = getPos(seg{1}, seg{2}, seg{3}, seg{4}, seg{5}, th_v(k), mu);
        
        % Aggiorna satellite
        set(h_sat, 'XData', r(1), 'YData', r(2), 'ZData', r(3));
        
        % Aggiorna traccia corrente
        t_x(end+1) = r(1); t_y(end+1) = r(2); t_z(end+1) = r(3);
        set(h_trail, 'XData', t_x, 'YData', t_y, 'ZData', t_z);
        
        drawnow;
        pause(0.005); % Regola velocità
    end
    
    % Marker Manovra
    plot3(t_x(end), t_y(end), t_z(end), 'x', 'MarkerEdgeColor', seg{9}, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
end
legend('show', 'Location', 'bestoutside', 'FontSize', 9);

