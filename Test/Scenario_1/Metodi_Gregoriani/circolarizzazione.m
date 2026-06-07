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

    T  = load("DatiSC1-2026.txt");
    gruppo = 23;
    
    % --- Orbita target / finale ---
    rr_f   = T(gruppo, 8:10)';
    vv_f   = T(gruppo, 11:13)';
    mu = 398600; 
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
% 2. OTTIMIZZAZIONE ORBITA AUSILIARIA CIRCOLARE (e_aux = 0)
% -------------------------------------------------------------------------
r_p_i = a_i*(1-e_i);                    % raggio pericentro orbita iniziale
r_p_f = a_f*(1-e_f);                    % raggio pericentro orbita finale

% L'orbita ausiliaria deve essere più grande dei pericentri per avere senso
% nei trasferimenti bitangenti 'pa' e 'ap'.
a_aux_min = max(r_p_i, r_p_f) + 500;    % a_aux_min = (r_p_i + r_a_f)/2 nel codice buttato;
a_aux_max = 900000;                     
N_search  = 5000;
a_aux_vec = linspace(a_aux_min, a_aux_max, N_search);
dv_total_vec = zeros(1, N_search);
dt_total_vec = zeros(1, N_search);

dv_salita_vec=zeros(2,N_search);
dv_discesa_vec=zeros(2,N_search);

% --- contributi separati al DeltaV ---
dv_plane_vec    = zeros(1, N_search);
dv_transfer_vec = zeros(1, N_search);

disp('Scansione di a_aux (circolare) in corso...');

% Vogliamo orbita finale circolare
e_aux_k = 0; % FORZIAMO L'ORBITA CIRCOLARE

for k = 1:N_search
    a_aux_k = a_aux_vec(k);

    % coasting fino al pericentro su orbita iniziale
    dt_coast_1 = TOF(a_i, e_i, th_i, 0, mu);
    
    % Manovra 1: bitangente 'pa' (da pericentro iniziale ad apocentro/orbita circolare)
    [dv1_1, dv1_2, dt_bitang_1] = bitangentTransfer(a_i, e_i, a_aux_k, e_aux_k, 'pa', mu);

    % Il coasting è scritto dopo il cambio di piano
    
    % Manovra 2: cambio di piano
    [dv2, om_tmp, th_plane_base] = changeOrbitalPlane(a_aux_k, e_aux_k, i_i, OM_i, om_i, i_f, OM_f, mu);

% --- FIX: SELEZIONE DEL NODO FUTURO (Rispetto a M2 che è a theta = pi) ---
% Calcoliamo i due nodi assicurandoci che siano nel range [0, 2*pi)
th_nodo1 = mod(th_plane_base, 2*pi);
th_nodo2 = mod(th_plane_base + pi, 2*pi);

% Il satellite parte da pi e va verso 2*pi. 
% Scegliamo il nodo che si trova in quel tratto di strada!
if th_nodo1 >= pi
    th_plane_tmp = th_nodo1;
else
    th_plane_tmp = th_nodo2;
end

    % A livello logico questo coasting è prima del cambio piano, ma ho
    % bisogno di calcolare th_plane_tmp prima
    % coasting da apocentro orbita circolare 1 fino a punto di cambio piano
    dt_coast_2 = TOF(a_aux_k, e_aux_k, pi, th_plane_tmp, mu);

    % Manovra 3: Rotazione pericentro: questa operazione deve dare un dv3
    % NULLO, perchè un'orbita circolare non ha motivo di essere ruotata: il
    % suo pericentro o apocentro è qualunque punto della circonferenza
    [dv3, th_3_i, th_3_f] = changePericenterArg(a_aux_k, e_aux_k, om_tmp, om_f, mu);
    
    % coasting per raggiungere il punto diametralmente opposto al pericentro finale
    dt_coast_3 = TOF(a_aux_k, e_aux_k, min(th_3_f), pi, mu);

    % Manovra 4: bitangente 'ap' (da apocentro/orbita circolare a pericentro finale)
    [dv4_1, dv4_2, dt_bitang_2] = bitangentTransfer(a_aux_k, e_aux_k, a_f, e_f, 'ap', mu);
    
% --- contributi separati ---
dv_transfer = abs(dv1_1) + abs(dv1_2) + abs(dv4_1) + abs(dv4_2);

dv_plane = abs(dv2);

% salvataggio vettori
dv_transfer_vec(k) = dv_transfer;
dv_plane_vec(k)    = dv_plane;

dv_salita_vec(:,k) = [dv1_1;dv1_2];
dv_discesa_vec(:,k)=[dv4_1;dv4_2];

% COSTO E TEMPO TOTALE
dv_total_vec(k) = dv_transfer + dv_plane + abs(dv3);

dt_total_vec(k) = dt_coast_1 + dt_bitang_1 + dt_coast_2 + dt_coast_3 + dt_bitang_2;
end

% Trova il minimo
[dv_min, idx_opt] = min(dv_total_vec);
a_aux = a_aux_vec(idx_opt);
e_aux = 0; % Circolare
dt = dt_total_vec(idx_opt);

fprintf('\n=== OTTIMIZZAZIONE COMPLETATA ===\n');
fprintf('  a_aux ottimale = %.4f km\n', a_aux);
fprintf('  e_aux          = %.6f (Circolare)\n', e_aux);
fprintf('  DeltaV TOTALE  = %.6f km/s\n', dv_min);
fprintf('  Deltat TOTALE  = %.6f gg\n', dt/86400);

% -------------------------------------------------------------------------
% PLOT ANDAMENTO DELTA-V
% -------------------------------------------------------------------------
figure('Name','Analisi DeltaV');

plot(a_aux_vec, dv_total_vec, 'LineWidth',2);
hold on;

plot(a_aux_vec, dv_transfer_vec, '--', 'LineWidth',1.8);

plot(a_aux_vec, dv_plane_vec, ':', 'LineWidth',2);

xline(a_aux,'r--','LineWidth',1.5);

grid on;
xlim([a_aux_min a_aux_max]);

xlabel('a_{aux} [km]');
ylabel('\DeltaV [km/s]');

title('Andamento del DeltaV');

legend('\DeltaV totale', ...
       '\DeltaV trasferimenti', ...
       '\DeltaV cambio piano', ...
       'a_{aux} ottimo', ...
       'Location','best');

figure;


plot(a_aux_vec,dv_salita_vec(1,:),a_aux_vec,dv_salita_vec(2,:), ...
    a_aux_vec,dv_discesa_vec(1,:),a_aux_vec,dv_discesa_vec(2,:), ...
    a_aux_vec,dv_transfer_vec, 'LineWidth',2);
legend('Salita 1','Salita 2','Discesa 1','Discesa 2','DV trasferimenti');
grid on;
title('DeltaV Traiettorie di Salita e Discesa');

% -------------------------------------------------------------------------
% 3. PLOT 3D COMPLETO
% -------------------------------------------------------------------------
figure('Name','Trasferimento Orbitale - Soluzione Circolare');

ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal; view(35, 25);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Sequenza di Manovre e Nodi (Orbita Ausiliaria Circolare)','FontSize',14);

% Terra con Mappa Topografica 
load topo; 
[XE, YE, ZE] = sphere(50);
surf(XE*6371, YE*6371, ZE*6371, 'CData', topo, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'HandleVisibility', 'off');
colormap(ax, topomap1); 

dth = deg2rad(1);

% --- PLOT ORBITE ---
% 1. Orbita Iniziale
plotOrbit_circolarizzazione(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [0.2 0.6 1], '-', 1.5, '1. Orbita Iniziale');
% 2. Trasferimento Salita
r_a_t1 = a_aux; 
a_t1 = (r_p_i + r_a_t1)/2;
e_t1 = (r_a_t1 - r_p_i)/(r_a_t1 + r_p_i);
plotOrbit_circolarizzazione(a_t1, e_t1, i_i, OM_i, om_i, 0, pi, dth, mu, [1 0.6 0.1], '-', 2, '2. Trasferimento Salita');
% 3. Orbita Ausiliaria (pre-piano)
plotOrbit_circolarizzazione(a_aux, e_aux, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [1 1 0], '--', 1.5, '3. Orbita Aux (Pre-Piano)');
% 4. Orbita Ausiliaria (post-piano)
plotOrbit_circolarizzazione(a_aux, e_aux, i_f, OM_f, om_tmp, 0, 2*pi, dth, mu, [0.8 0.2 0.8], '-.', 1.5, '4. Orbita Aux (Post-Piano)');
% 5. Trasferimento Discesa
a_t2 = (a_aux + r_p_f)/2;
e_t2 = (a_aux - r_p_f)/(a_aux + r_p_f);
plotOrbit_circolarizzazione(a_t2, e_t2, i_f, OM_f, om_f, pi, 2*pi, dth, mu, [1 0.3 0.3], ':', 2, '5. Trasferimento Discesa');
% 6. Orbita Finale Target
plotOrbit_circolarizzazione(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, dth, mu, [0 1 0], '-', 2, '6. Orbita Target');

% --- PLOT MARKER DELLE MANOVRE (Calcolati analiticamente) ---
% M1: Partenza dal pericentro iniziale (th = 0)
P1 = getPos3D(a_i, e_i, i_i, OM_i, om_i, 0); 
% M2: Arrivo in quota, apocentro ellisse di trasferimento (th = pi sull'orbita di parcheggio)
P2 = getPos3D(a_aux, e_aux, i_i, OM_i, om_i, pi); 
% M3: Cambio Piano, calcolato dalla funzione changeOrbitalPlane
P3 = getPos3D(a_aux, e_aux, i_i, OM_i, om_i, th_plane_tmp); 
% M4: Inizio discesa, diametralmente opposto al pericentro finale (th = pi nel sistema finale)
P4 = getPos3D(a_aux, e_aux, i_f, OM_f, om_f, pi); 
% M5: Arrivo al pericentro target (th = 0 nel sistema finale)
P5 = getPos3D(a_f, e_f, i_f, OM_f, om_f, 0); 

% Disegno i Marker
plot3(P1(1), P1(2), P1(3), 'wo', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'M1: Inizio Salita');
plot3(P2(1), P2(2), P2(3), 'wo', 'MarkerFaceColor', 'y', 'MarkerSize', 7, 'DisplayName', 'M2: Fine Salita');
plot3(P3(1), P3(2), P3(3), 'wo', 'MarkerFaceColor', 'm', 'MarkerSize', 7, 'DisplayName', 'M3: Cambio Piano');
plot3(P4(1), P4(2), P4(3), 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'M4: Inizio Discesa');
plot3(P5(1), P5(2), P5(3), 'wo', 'MarkerFaceColor', 'c', 'MarkerSize', 7, 'DisplayName', 'M5: Arrivo Target');

legend('show','Location','eastoutside','TextColor','w','Color','k','FontSize',11);
hold off;

% =========================================================================
% DA INCOLLARE ALLA FINE - RICERCA MINIMO LOCALE CON FMINBND (CORRETTO)
% =========================================================================

% 1. Definizione dei parametri geometrici
r1 = r_p_i;
r2 = r_p_f;

% Calcolo dell'angolo di cambio piano effettivo (alpha)
cos_alpha = cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(OM_f - OM_i);
alpha = acos(max(-1, min(1, cos_alpha))); 

% 2. Definizione della funzione di costo CON I VALORI ASSOLUTI (Cruciale!)
dv_fun = @(x) ...
    abs(sqrt(2*mu/r1 - 2*mu./(r1+x)) - sqrt(2*mu/r1 - mu/a_i)) ... % M1: Salita
  + abs(sqrt(mu./x) - sqrt(2*mu./x - 2*mu./(r1+x))) ...            % M2: Parcheggio 1
  + abs(2 * sqrt(mu./x) * sin(alpha/2)) ...                        % M3: Cambio piano
  + abs(sqrt(mu./x) - sqrt(2*mu./x - 2*mu./(r2+x))) ...            % M4: Inizio discesa
  + abs(sqrt(2*mu/r2 - 2*mu./(r2+x)) - sqrt(2*mu/r2 - mu/a_f));    % M5: Arrivo

% 3. Esecuzione fminbnd con confini stretti
% Cerchiamo solo nella prima porzione, prima della discesa asintotica
limite_inf = a_aux_min;
limite_sup = 100000; 

% Opzioni silenziose (Display = off per non stampare a schermo)
opzioni = optimset('Display', 'off', 'TolX', 1e-8);

[a_aux_ottimo, dv_minimo_analitico] = fminbnd(dv_fun, limite_inf, limite_sup, opzioni);

fprintf('\n========================================================\n');
fprintf('  OTTIMIZZAZIONE FMINBND (Risultato Analitico)\n');
fprintf('========================================================\n');
fprintf('  a_aux locale            : %.4f km\n', a_aux_ottimo);
fprintf('  DeltaV minimo analitico : %.6f km/s\n', dv_minimo_analitico);
fprintf('========================================================\n');

% 4. Aggiunta del marker al plot esistente
fig_dv = findobj('type', 'figure', 'Name', 'Analisi DeltaV');
if ~isempty(fig_dv)
    figure(fig_dv); hold on;
    plot(a_aux_ottimo, dv_minimo_analitico, 'p', 'MarkerSize', 14, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 1.5, ...
        'DisplayName', 'Minimo Locale (Esatto)');
    legend('Location', 'best');
end
%% =========================================================================
% 4. ANIMAZIONE 3D AVANZATA - TRACCIA MULTI-COLORE E ORBITE DI RIFERIMENTO
% =========================================================================

% --- 1. SETUP SCENA ---
fig_anim = figure('Name', 'Simulazione Dinamica Trasferimento Orbitale Circolare', ...
    'Units','normalized','Position',[0.1 0.1 0.8 0.8]);

% Assi scuri ad alto contrasto
ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal; view(35, 25);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Animazione Sequenza Manovre (Orbita Ausiliaria Circolare)', 'Color','w','FontSize', 14);

% Terra 3D
[xE, yE, zE] = sphere(50);
try
    load topo topo topomap1;
    surf(xE * 6371, yE * 6371, zE * 6371, 'FaceColor', 'texturemap', 'CData', topo, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    colormap(ax, topomap1);
catch
    surf(xE * 6371, yE * 6371, zE * 6371, 'FaceColor', [0.1 0.4 0.8], 'EdgeColor', 'none');
end

% --- 2. PLOT ORBITE DI RIFERIMENTO (Tratteggiate per la legenda) ---
plotOrbit_circolarizzazioneStaticAnim(a_i, e_i, i_i, OM_i, om_i, mu, [0.2 0.6 1], '--', 1, 'Rif. Iniziale');
plotOrbit_circolarizzazioneStaticAnim(a_aux, e_aux, i_i, OM_i, om_i, mu, [1 1 0], '--', 1, 'Rif. Parcheggio Pre-Piano');
plotOrbit_circolarizzazioneStaticAnim(a_aux, e_aux, i_f, OM_f, om_f, mu, [0.8 0.2 0.8], '--', 1, 'Rif. Parcheggio Post-Piano');
plotOrbit_circolarizzazioneStaticAnim(a_f, e_f, i_f, OM_f, om_f, mu, [0 1 0], '--', 1, 'Rif. Finale Target');

% --- 3. DEFINIZIONE SEGMENTI ---
num_punti_curva = 500; 

segmenti = {};
% 1. Coasting iniziale verso il pericentro
th_s = th_i; th_e = 0; if th_e <= th_s, th_e = th_e + 2*pi; end
segmenti{1} = {a_i, e_i, i_i, OM_i, om_i, th_s, th_e, 'Fase 1: Approccio', [0.2 0.6 1]};

% 2. Manovra di Salita
segmenti{2} = {a_t1, e_t1, i_i, OM_i, om_i, 0, pi, 'Fase 2: Manovra Salita', [1 0.6 0.1]};

% 3. Coasting su Aux (Parcheggio Pre-Piano)
th_s = pi; th_e = th_plane_tmp; if th_e <= th_s, th_e = th_e + 2*pi; end
segmenti{3} = {a_aux, e_aux, i_i, OM_i, om_i, th_s, th_e, 'Fase 3: Parcheggio Pre-Piano', [1 1 0]};

% 4. Allineamento dopo cambio piano
th_post_plane_new = mod(th_plane_tmp + om_tmp - om_f, 2*pi);
th_s = th_post_plane_new; th_e = pi; if th_e <= th_s, th_e = th_e + 2*pi; end
segmenti{4} = {a_aux, e_aux, i_f, OM_f, om_f, th_s, th_e, 'Fase 4: Parcheggio Post-Piano', [0.8 0.2 0.8]};

% 5. Manovra di Discesa
segmenti{5} = {a_t2, e_t2, i_f, OM_f, om_f, pi, 2*pi, 'Fase 5: Manovra Discesa', [1 0.3 0.3]};

% 6. Arrivo a Target
th_s = 0; th_e = th_f; if th_e <= th_s, th_e = th_e + 2*pi; end
segmenti{6} = {a_f, e_f, i_f, OM_f, om_f, th_s, th_e, 'Fase 6: Inserimento Target', [0 1 0]};

% --- 4. CICLO DI ANIMAZIONE ---
h_sat = plot3(NaN, NaN, NaN, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Satellite');

% Legenda ad alto contrasto
legend('show', 'Location', 'bestoutside', 'FontSize', 9, 'TextColor', 'w', 'Color', 'k', 'EdgeColor', [0.5 0.5 0.5]);

% === CONTROLLO VELOCITÀ (LA VIA DI MEZZO) ===
skip_frames = 5; 

for s = 1:length(segmenti)
    seg = segmenti{s};
    th_v = linspace(seg{6}, seg{7}, num_punti_curva); 
    
    h_trail = plot3(NaN, NaN, NaN, 'Color', seg{9}, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    t_x = NaN(1, num_punti_curva); 
    t_y = NaN(1, num_punti_curva); 
    t_z = NaN(1, num_punti_curva);
    
    for k = 1:num_punti_curva
        if ~isvalid(h_sat), return; end 
        
        r = getPosAnim(seg{1}, seg{2}, seg{3}, seg{4}, seg{5}, th_v(k), mu);
        
        t_x(k) = r(1); 
        t_y(k) = r(2); 
        t_z(k) = r(3);
        
        if mod(k, skip_frames) == 0 || k == num_punti_curva
            set(h_sat, 'XData', r(1), 'YData', r(2), 'ZData', r(3));
            set(h_trail, 'XData', t_x, 'YData', t_y, 'ZData', t_z);
            drawnow; 
            pause(0.001); % Freno a mano per stabilizzare la fluidità
        end
    end
    % Marker delle manovre (nascosti dalla legenda)
    plot3(t_x(k), t_y(k), t_z(k), 'x', 'MarkerEdgeColor', seg{9}, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
end

