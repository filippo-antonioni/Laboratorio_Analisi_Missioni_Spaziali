% =========================================================================
% MASTER SCRIPT - SCENARIO 1 (Ricerca Globale + TOF + Grafica Pulita)
% Usa SOLO le funzioni utente originali
% =========================================================================
clear; clc; close all;
disp('--- AVVIO RICERCA GLOBALE CON DOPPIO CICLO FOR ---');

% --- 1. CARICAMENTO DATI ---
try
    T = load("DatiSC1-2026.txt");
    gruppo = 23;
    rr_f = T(gruppo,8:10)'; vv_f = T(gruppo,11:13)';
    orbita_iv = T(gruppo,2:7);
catch
    rr_f = [-11441; -7209.8; -1302.9]; vv_f = [1.214; -1.711; -4.716];
    orbita_iv = [24400, 0.7283, 0.1047, 2.894, 3.107, 2.011];
end

mu = 398600;
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);
a_i = orbita_iv(1); e_i = orbita_iv(2); i_i = orbita_iv(3);
OM_i = orbita_iv(4); om_i = orbita_iv(5); th_i = orbita_iv(6);

% --- 2. SETUP DEL DOPPIO CICLO FOR ---
r_p_i = a_i * (1 - e_i); 
r_a_i = a_i * (1 + e_i); 
r_p_f = a_f * (1 - e_f); 

r_a_max = 350000; 
r_a_test_vec = linspace(max(r_a_i, r_p_f), r_a_max, 100); 

min_costo_A = inf; 
min_costo_B = inf; 

fprintf('Esplorazione griglia (Apocentri e Pericentri) in corso...\n');

for idx_a = 1:length(r_a_test_vec)
    r_a_park = r_a_test_vec(idx_a);
    r_p_min = 6371 + 300; 
    r_p_test_vec = linspace(r_p_min, r_a_park - 1000, 100); 
    
    for idx_p = 1:length(r_p_test_vec)
        r_p_park = r_p_test_vec(idx_p);
        
        a_park = (r_a_park + r_p_park) / 2;
        e_park = (r_a_park - r_p_park) / (r_a_park + r_p_park);
        if e_park >= 1 || e_park < 0, continue; end
        
        % -------------------------------------------------------------
        % FASI COMUNI 
        % -------------------------------------------------------------
        % Coasting iniziale fino al pericentro (th=0)
        dt_coast1 = TOF(a_i, e_i, th_i, 0, mu);
        
        % Bitangente 1 (Salita 'pa')
        [dV1_A, dV1_B, dt1] = bitangentTransfer(a_i, e_i, a_park, e_park, 'pa', mu);
        costo_salita = abs(dV1_A) + abs(dV1_B);
        
        % Cambio Piano 
        [dV_plane, om_fun, theta_plane] = changeOrbitalPlane(a_park, e_park, i_i, OM_i, om_i, i_f, OM_f, mu);
        costo_piano = abs(dV_plane);
        
        % Attesa verso il nodo del cambio piano
        dt_coast_park1 = TOF(a_park, e_park, pi, theta_plane, mu);
        
        % -------------------------------------------------------------
        % SEQUENZA A: Discesa ('ap') -> Cambio Pericentro
        % -------------------------------------------------------------
        [dV_arg_A, thi_fun_A, thf_fun_A] = changePericenterArg(a_f, e_f, om_fun, om_f, mu);
        thi_A = thi_fun_A(1); thf_A = thf_fun_A(1); 
        
        [dV2A_A, dV2A_B, dt2A] = bitangentTransfer(a_park, e_park, a_f, e_f, 'ap', mu);
        
        dt_coast_park2_A = TOF(a_park, e_park, theta_plane, pi, mu); 
        dt_coast_fin1_A = TOF(a_f, e_f, 0, thi_A, mu); 
        dt_coast_fin2_A = TOF(a_f, e_f, thf_A, th_f, mu); 
        
        tof_A = dt_coast1 + dt1 + dt_coast_park1 + dt_coast_park2_A + dt2A + dt_coast_fin1_A + dt_coast_fin2_A;
        totale_A = costo_salita + costo_piano + abs(dV2A_A) + abs(dV2A_B) + abs(dV_arg_A);
        
        if totale_A < min_costo_A
            min_costo_A = totale_A;
            best_A = struct('a_park', a_park, 'e_park', e_park, 'costo', totale_A, 'tof', tof_A);
        end
        
        % -------------------------------------------------------------
        % SEQUENZA B: Cambio Pericentro -> Discesa ('ap')
        % -------------------------------------------------------------
        [dV_arg_B, thi_fun_B, thf_fun_B] = changePericenterArg(a_park, e_park, om_fun, om_f, mu);
        thi_B = thi_fun_B(1); thf_B = thf_fun_B(1);
        
        [dV2B_A, dV2B_B, dt2B] = bitangentTransfer(a_park, e_park, a_f, e_f, 'ap', mu);
        
        dt_coast_park2_B = TOF(a_park, e_park, theta_plane, thi_B, mu); 
        dt_coast_park3_B = TOF(a_park, e_park, thf_B, pi, mu); 
        dt_coast_fin_B = TOF(a_f, e_f, 0, th_f, mu); 
        
        tof_B = dt_coast1 + dt1 + dt_coast_park1 + dt_coast_park2_B + dt_coast_park3_B + dt2B + dt_coast_fin_B;
        totale_B = costo_salita + costo_piano + abs(dV_arg_B) + abs(dV2B_A) + abs(dV2B_B);
        
        if totale_B < min_costo_B
            min_costo_B = totale_B;
            best_B = struct('a_park', a_park, 'e_park', e_park, 'costo', totale_B, 'tof', tof_B, ...
                            'dv1', dV1_A, 'dv2', dV1_B, 'dv_p', dV_plane, ...
                            'dv_a', dV_arg_B, 'dv3', dV2B_A, 'dv4', dV2B_B, ...
                            'om_fun', om_fun, 'th_plane', theta_plane, 'th_peri', thi_B, ...
                            't_coast1', dt_coast1, 't_trsf1', dt1, 't_att1', dt_coast_park1, ...
                            't_att2', dt_coast_park2_B, 't_att3', dt_coast_park3_B, ...
                            't_trsf2', dt2B, 't_coastF', dt_coast_fin_B);
        end
    end
end

% --- 3. STAMPA RISULTATI E CONFRONTO ---
fprintf('\n=======================================================\n');
fprintf('CONFRONTO SEQUENZE (Costi e TOF)\n');
fprintf('=======================================================\n');
fprintf('Sequenza A (Discesa -> Pericentro): %.4f km/s | TOF: %.2f giorni\n', min_costo_A, best_A.tof/86400);
fprintf('Sequenza B (Pericentro -> Discesa): %.4f km/s | TOF: %.2f giorni\n', min_costo_B, best_B.tof/86400);

if min_costo_B < min_costo_A
    fprintf('\n>>> VITTORIA SEQUENZA B! Invertire le manovre fa risparmiare %.4f km/s!\n', min_costo_A - min_costo_B);
    opt = best_B; 
else
    fprintf('\n>>> VITTORIA SEQUENZA A!\n');
    opt = best_A;
end

fprintf('\n--- DETTAGLIO MANOVRE (Sequenza Vincitrice) ---\n');
fprintf('Raggio Apocentro Appoggio:  %.2f km\n', opt.a_park*(1+opt.e_park));
fprintf('Raggio Pericentro Appoggio: %.2f km\n', opt.a_park*(1-opt.e_park));
fprintf('Impulso 1 (Inizio Salita):  %.4f km/s\n', abs(opt.dv1));
fprintf('Impulso 2 (Fine Salita):    %.4f km/s\n', abs(opt.dv2));
fprintf('Impulso 3 (Cambio Piano):   %.4f km/s\n', abs(opt.dv_p));
fprintf('Impulso 4 (Cambio Peri):    %.4f km/s\n', abs(opt.dv_a));
fprintf('Impulso 5 (Inizio Discesa): %.4f km/s\n', abs(opt.dv3));
fprintf('Impulso 6 (Fine Discesa):   %.4f km/s\n', abs(opt.dv4));
fprintf('TOTALE ASSOLUTO:            %.4f km/s\n', opt.costo);

fprintf('\n--- TEMPI DI VOLO (TOF) ---\n');
fprintf('Coasting a Pericentro Iniziale:   %.2f giorni\n', opt.t_coast1 / 86400);
fprintf('Trasferimento 1 (Salita):         %.2f giorni\n', opt.t_trsf1 / 86400);
fprintf('Coasting a Cambio Piano:          %.2f giorni\n', opt.t_att1 / 86400);
fprintf('Coasting a Cambio Pericentro:     %.2f giorni\n', opt.t_att2 / 86400);
fprintf('Coasting a Discesa (Apocentro):   %.2f giorni\n', opt.t_att3 / 86400);
fprintf('Trasferimento 2 (Discesa):        %.2f giorni\n', opt.t_trsf2 / 86400);
fprintf('Coasting a Posizione Finale:      %.2f giorni\n', opt.t_coastF / 86400);
fprintf('TEMPO TOTALE DI MISSIONE:         %.2f giorni\n', opt.tof / 86400);

% =========================================================================
% --- 4. GRAFICA AVANZATA 3D (Finestra Classica) ---
% =========================================================================
fig = figure('Name', 'Analisi Missione: Trasferimento Ottimizzato', 'NumberTitle', 'off');
hold on; grid on; axis equal; view(3); rotate3d on;
xlabel('X [km]', 'FontWeight', 'bold'); 
ylabel('Y [km]', 'FontWeight', 'bold'); 
zlabel('Z [km]', 'FontWeight', 'bold');
title('Traiettoria 3D Ottimizzata', 'FontSize', 14, 'FontWeight', 'bold');

th_vec = linspace(0, 2*pi, 300);
get_c = @(a,e,i,OM,om,th,idx) subsref(par2car(a,e,i,OM,om,th,mu), struct('type','()','subs',{{idx}}));

r_p1 = a_i*(1-e_i); r_a1 = opt.a_park*(1+opt.e_park);
a_t1 = (r_p1 + r_a1)/2; e_t1 = (r_a1 - r_p1)/(r_a1 + r_p1);
r_a2 = opt.a_park*(1+opt.e_park); r_p2 = a_f*(1-e_f);
a_t2 = (r_a2 + r_p2)/2; e_t2 = (r_a2 - r_p2)/(r_a2 + r_p2);

% --- DISEGNO ORBITE ---
plot3(arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,1), th_vec), arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,2), th_vec), arrayfun(@(th) get_c(a_i,e_i,i_i,OM_i,om_i,th,3), th_vec), 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', 'Orbita Iniziale');
plot3(arrayfun(@(th) get_c(a_t1,e_t1,i_i,OM_i,om_i,th,1), th_vec), arrayfun(@(th) get_c(a_t1,e_t1,i_i,OM_i,om_i,th,2), th_vec), arrayfun(@(th) get_c(a_t1,e_t1,i_i,OM_i,om_i,th,3), th_vec), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'Trsf 1: Salita');
plot3(arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_i,OM_i,om_i,th,1), th_vec), arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_i,OM_i,om_i,th,2), th_vec), arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_i,OM_i,om_i,th,3), th_vec), 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Pre-Plane)');
plot3(arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_f,OM_f,opt.om_fun,th,1), th_vec), arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_f,OM_f,opt.om_fun,th,2), th_vec), arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_f,OM_f,opt.om_fun,th,3), th_vec), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Post-Plane)');
plot3(arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_f,OM_f,om_f,th,1), th_vec), arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_f,OM_f,om_f,th,2), th_vec), arrayfun(@(th) get_c(opt.a_park,opt.e_park,i_f,OM_f,om_f,th,3), th_vec), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Post-Peri)');
plot3(arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,om_f,th,1), th_vec), arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,om_f,th,2), th_vec), arrayfun(@(th) get_c(a_t2,e_t2,i_f,OM_f,om_f,th,3), th_vec), 'Color', [0.6350 0.0780 0.1840], 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'Trsf 2: Discesa');
plot3(arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,1), th_vec), arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,2), th_vec), arrayfun(@(th) get_c(a_f,e_f,i_f,OM_f,om_f,th,3), th_vec), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.5, 'DisplayName', 'Orbita Finale');

% --- MARKERS DELLE MANOVRE ---
[r1,~] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu);
plot3(r1(1), r1(2), r1(3), 'ok', 'MarkerSize', 7, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'DisplayName', '1. Inizio Salita');
[r2,~] = par2car(a_t1, e_t1, i_i, OM_i, om_i, pi, mu);
plot3(r2(1), r2(2), r2(3), 'sk', 'MarkerSize', 7, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'DisplayName', '2. Arrivo su Appoggio');
[r3,~] = par2car(opt.a_park, opt.e_park, i_i, OM_i, om_i, opt.th_plane, mu);
plot3(r3(1), r3(2), r3(3), '^k', 'MarkerSize', 8, 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', '3. Cambio Piano');
[r4,~] = par2car(opt.a_park, opt.e_park, i_f, OM_f, opt.om_fun, opt.th_peri, mu);
plot3(r4(1), r4(2), r4(3), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'DisplayName', '4. Cambio Pericentro');
[r5,~] = par2car(opt.a_park, opt.e_park, i_f, OM_f, om_f, pi, mu);
plot3(r5(1), r5(2), r5(3), 'dk', 'MarkerSize', 7, 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'DisplayName', '5. Inizio Discesa');
[r6,~] = par2car(a_t2, e_t2, i_f, OM_f, om_f, 0, mu);
plot3(r6(1), r6(2), r6(3), 'hk', 'MarkerSize', 8, 'MarkerFaceColor', [0.4660 0.6740 0.1880], 'DisplayName', '6. Arrivo Orbita Finale');

% --- TERRA 3D FOTOREALISTICA ---
R_earth = 6371;
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

% =========================================================================
% FUNZIONI ORIGINALI (Esattamente come fornite dall'utente)
% =========================================================================
