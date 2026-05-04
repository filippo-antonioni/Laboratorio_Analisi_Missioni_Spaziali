
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

% --- 2. SETUP DEL DOPPIO CICLO FOR ---
r_p_i = a_i * (1 - e_i); 
r_a_i = a_i * (1 + e_i); 
r_p_f = a_f * (1 - e_f); 

% INIZIO DELLA SEQUENZA DI MANOVRE

%coasting dal punto di partenza fino all'apocentro
dt_coast1 = TOF(a_i, e_i, th_i, pi, mu);

%calcolo posizioni e velocità all'apocentro sull'orbita iniziale 
[rr_i_a, vv_i_a] = par2car(a_i, e_i, i_i, OM_i, om_i, pi, mu);
%calcolo posizione e velocità all'apocentro nell'orbita di parcheggio
%circolare
[rr_park, vv_park] = par2car(a_i, 0, i_i, OM_i, om_i, pi, mu);

costo_circolarizzazione = norm(vv_i_a-vv_park);
a_park=a_i;
e_park=0;

%calcolo deltav del cambio piano
[dV_plane, om_plane, theta_plane] = changeOrbitalPlane(a_park, e_park, i_i, OM_i, om_i, i_f, OM_f, mu);
costo_piano = abs(dV_plane);
        
 th_plane_tmp = theta_plane;
 %coasting per arrivare al punto di cambio piano
 dt_coast_park = TOF(a_park, e_park, pi, th_plane_tmp, mu);
        
 % Cambio anomalia pericentro sull'orbita di appoggio
 [dV_arg, thi_fun, thf_fun] = changePericenterArg(a_park, e_park, om_plane, om_f, mu);
 
 %sono abbastanza convinto del fatto che fare un coastin per arrivare al
 %cambio del pericentro abbia poco senso visto che non è una vera manovra
 %in questo caso b
 [thi_arg_max, idx_thi_arg_max] = max(thi_fun);
 [thi_arg_min, idx_thi_arg_min] = min(thi_fun);
        
        % Calcolo tempo coasting da punto cambio piano al punto cambio pericentro
        if theta_plane > thi_arg_min && theta_plane < thi_arg_max
            dt_coast_periarg = TOF(a_park, e_park, theta_plane, thi_arg_max, mu);
            arg = 1;
            thi_arg = thi_arg_max;
        else 
            dt_coast_periarg = TOF(a_park, e_park, theta_plane, thi_arg_min, mu); 
            arg = 0;
            thi_arg = thi_arg_min;
        end
        
        switch arg
            case 1
                thf_arg=thf_fun(idx_thi_arg_max);
            
            otherwise
                thf_arg=thf_fun(idx_thi_arg_min);
            
        end

        % Calcolo coasting dal punto post-cambio anomalia fino all'apocentro
        % if arg
        %     dt_coast_park2 = TOF(a_park, e_park, thf_fun(idx_thi_arg_max), pi, mu); 
        % else 
        %     dt_coast_park2 = TOF(a_park, e_park, thf_fun(idx_thi_arg_min), pi, mu); 
        % end
       
        % provo diverse alternative 
        
        % OPZIONE 1 DISCESA CON BITANGENTE
                
        
        dv_bitang=[];
        
        [DeltaV1_ap, DeltaV2_ap, Deltat_ap] = bitangentTransfer(a_park, e_park, a_f, e_f, 'ap', mu);
        dv_bitang=[dv_bitang; abs(DeltaV1_ap)+abs(DeltaV2_ap)];

        [DeltaV1_pa, DeltaV2_pa, Deltat_pa] = bitangentTransfer(a_park, e_park, a_f, e_f, 'pa', mu);
        dv_bitang=[dv_bitang; abs(DeltaV1_pa)+abs(DeltaV2_pa)];

        [DeltaV1_aa, DeltaV2_aa, Deltat_aa] = bitangentTransfer(a_park, e_park, a_f, e_f, 'aa', mu);
        dv_bitang=[dv_bitang; abs(DeltaV1_aa)+abs(DeltaV2_aa)];

        [DeltaV1_pp, DeltaV2_pp, Deltat_pp] = bitangentTransfer(a_park, e_park, a_f, e_f, 'pp', mu);
        dv_bitang=[dv_bitang; abs(DeltaV1_pp)+abs(DeltaV2_pp)];
        
        [dv_bitang2,idx_bitang2]=min(dv_bitang);
        
        switch idx_bitang2
            case 1
                th_bitang=pi;
                th_arrivo_bitang=0;
                dt_bitang2=Deltat_ap;
            case 2
                th_bitang=0;
                th_arrivo_bitang=pi;
                dt_bitang2=Deltat_pa;
            case 3
                th_bitang=pi;
                th_arrivo_bitang=pi;
                dt_bitang2=Deltat_aa;
            case 4
                th_bitang=0;
                th_arrivo_bitang=0;
                dt_bitang2=Deltat_pp;
            otherwise 
                error('\nDimensione vettore non valida');
        end
        
        %coasting dal punto di cambio anomalia del pericentro che secondo
        %me va bene ovunque fino al punto della bitangente scelta 
        dt_coast_park2=TOF(a_park,e_park,thf_arg,th_bitang,mu);
        
        dt_coast_finale_A=TOF(a_f,e_f,th_arrivo_bitang,th_f,mu);  
          
        tempo_tot_A = dt_coast1 + dt_coast_park + dt_coast_periarg + dt_bitang2 + dt_coast_park2+ dt_coast_finale_A;

        dV_tot_A = costo_circolarizzazione + costo_piano + abs(dV_arg) + dv_bitang2;

        % OPZIONE 2 FACCIO MANOVRA DIRETTAMENTE DA ORBITA POST PLANE POST
        % ARG A QUELLA FINALE SFRUTTANDO INTERSEZIONE CHE E' CERTA

        [interseca, th1_int, th2_int] = orbit_intersection(a_park, e_park, om_f, a_f, e_f, om_f);
        if interseca 
            
            if thf_arg< th1_int(1) || thf_arg>th1_int(2)
                
                dt_coast_intersez=TOF(a_park,e_park,thf_arg,th1_int(1),mu);
                th1_intersez=th1_int(1);
                th2_intersez=th2_int(1);
            else 
                dt_coast_intesez=TOF(a_park,e_park,thf_arg,th1_int(2),mu);
                th1_intersez=th1_int(2);
                th2_intersez=th2_int(2);
            end
                
            [rr_intersez_1, vv_intersez_1] = par2car(a_park, e_park, i_f, OM_f, om_f, th1_intersez, mu);
            [rr_intersez_2, vv_intersez_2] = par2car(a_f, e_f, i_f, OM_f, om_f, th2_intersez, mu); 
            dv_discesa=norm(vv_intersez_1-vv_intersez_2);
            dt_coast_finale_B=TOF(a_f,e_f,th2_intersez,th_f,mu);

            tempo_tot_B = dt_coast1 + dt_coast_park + dt_coast_periarg + dt_coast_intersez + dt_coast_finale_B;

             dV_tot_B = costo_circolarizzazione + costo_piano + abs(dV_arg) + dv_discesa;
            
        else 
            error('\nOrbite non si intesecano'); 
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
fprintf('   CONFRONTO STRATEGIE: OPZIONE A vs OPZIONE B\n');
fprintf('=======================================================\n');

% --- STAMPA OPZIONE A ---
fprintf('\n---> OPZIONE A: Discesa con Bitangente\n');
fprintf('DeltaV Totale A: %.4f km/s | TOF: %.2f giorni\n', dV_tot_A, tempo_tot_A/86400);
fprintf('- Impulso 1 (Circolarizzazione):%.4f km/s\n', costo_circolarizzazione);
fprintf('- Impulso 2 (Cambio Piano):     %.4f km/s\n', costo_piano);
fprintf('- Impulso 3 (Cambio Pericentro):%.4f km/s\n', abs(dV_arg));
fprintf('- Impulsi 4+5 (Bitangente):     %.4f km/s\n', dv_bitang2);
fprintf('  [Tempi di Volo A]\n');
fprintf('   Coasting Iniziale:      %.2f gg\n', dt_coast1 / 86400);
fprintf('   Coasting Cambio Piano:  %.2f gg\n', dt_coast_park / 86400);
fprintf('   Coasting Cambio Peri:   %.2f gg\n', dt_coast_periarg / 86400);
fprintf('   Coasting a Bitangente:  %.2f gg\n', dt_coast_park2 / 86400);
fprintf('   Trasferimento Bitang:   %.2f gg\n', dt_bitang2 / 86400);
fprintf('   Coasting Finale:        %.2f gg\n', dt_coast_finale_A / 86400);

% --- STAMPA OPZIONE B ---
% NOTA: Calcolo il tempo di intersezione per differenza per evitare
% eventuali "typo" di battitura nelle variabili del tuo blocco if/else
dt_coast_intersez_print = tempo_tot_B - dt_coast1 - dt_coast_park - dt_coast_periarg - dt_coast_finale_B;

fprintf('\n---> OPZIONE B: Intersezione Diretta\n');
fprintf('DeltaV Totale B: %.4f km/s | TOF: %.2f giorni\n', dV_tot_B, tempo_tot_B/86400);
fprintf('- Impulso 1 (Circolarizzazione):%.4f km/s\n', costo_circolarizzazione);
fprintf('- Impulso 2 (Cambio Piano):     %.4f km/s\n', costo_piano);
fprintf('- Impulso 3 (Cambio Pericentro):%.4f km/s\n', abs(dV_arg));
fprintf('- Impulso 4 (Salto Diretto):    %.4f km/s\n', dv_discesa);
fprintf('  [Tempi di Volo B]\n');
fprintf('   Coasting Iniziale:      %.2f gg\n', dt_coast1 / 86400);
fprintf('   Coasting Cambio Piano:  %.2f gg\n', dt_coast_park / 86400);
fprintf('   Coasting Cambio Peri:   %.2f gg\n', dt_coast_periarg / 86400);
fprintf('   Coasting a Intersez:    %.2f gg\n', dt_coast_intersez_print / 86400);
fprintf('   Coasting Finale:        %.2f gg\n', dt_coast_finale_B / 86400);

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

% Circolarizzazione all'Apocentro (Nessun ellisse di trasferimento, sale diretta)
[r_circ,~] = par2car(a_i, e_i, i_i, OM_i, om_i, pi, mu);
plot3(r_circ(1), r_circ(2), r_circ(3), 'sk', 'MarkerSize', 7, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'DisplayName', 'Circolarizzazione');

% 2. Orbita Appoggio Pre-Piano
plot3(arrayfun(@(th) get_c(a_park,e_park,i_i,OM_i,om_i,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_i,OM_i,om_i,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_i,OM_i,om_i,th,3), th_vec_full), ...
      'Color', [0.9290 0.6940 0.1250], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Pre-Plane)');

[r_plane,~] = par2car(a_park, e_park, i_i, OM_i, om_i, theta_plane, mu);
plot3(r_plane(1), r_plane(2), r_plane(3), '^k', 'MarkerSize', 8, 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Cambio Piano');

% 3. Orbita Appoggio Post-Piano
plot3(arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_plane,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_plane,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_plane,th,3), th_vec_full), ...
      'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Post-Plane)');

[r_cambio_peri,~] = par2car(a_park, e_park, i_f, OM_f, om_plane, thi_arg, mu);
plot3(r_cambio_peri(1), r_cambio_peri(2), r_cambio_peri(3), 'pk', 'MarkerSize', 10, 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'DisplayName', 'Cambio Pericentro');

% 4. Orbita Appoggio Post-Peri
plot3(arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_f,th,1), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_f,th,2), th_vec_full), ...
      arrayfun(@(th) get_c(a_park,e_park,i_f,OM_f,om_f,th,3), th_vec_full), ...
      'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Appoggio (Post-Peri)');

% --- DISEGNO DINAMICO DELLA MANOVRA FINALE IN BASE A CHI HA VINTO ---
if opzione(1) == 'A'
    % Ricostruisco l'ellisse della bitangente per plottarla in 3D
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
else
    % Se ha vinto l'intersezione diretta metto solo il marker del salto
    if exist('th1_intersez', 'var')
        [r_intersez,~] = par2car(a_park, e_park, i_f, OM_f, om_f, th1_intersez, mu);
        plot3(r_intersez(1), r_intersez(2), r_intersez(3), 'hk', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Punto di Intersezione (B)');
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