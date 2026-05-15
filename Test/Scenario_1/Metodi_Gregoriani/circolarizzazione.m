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

% -------------------------------------------------------------------------
% 2. OTTIMIZZAZIONE ORBITA AUSILIARIA CIRCOLARE (e_aux = 0)
% -------------------------------------------------------------------------
r_p_i = a_i*(1-e_i);                    % raggio pericentro orbita iniziale
r_p_f = a_f*(1-e_f);                    % raggio pericentro orbita finale

% L'orbita ausiliaria deve essere più grande dei pericentri per avere senso
% nei trasferimenti bitangenti 'pa' e 'ap'.
a_aux_min = max(r_p_i, r_p_f) + 500;    % a_aux_min = (r_p_i + r_a_f)/2 nel codice buttato;
a_aux_max = 320000;                     
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

dv_salita_vec(:,k) = [abs(dv1_1);abs(dv1_2)];
dv_discesa_vec(:,k)=[abs(dv4_1);abs(dv4_2)];

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
    a_aux_vec,dv_transfer_vec);
legend('Salita 1','Salita 2','Discesa 1','Discesa 2','DV trasferimenti');
grid on;
title('DeltaV salita e discesa');

% -------------------------------------------------------------------------
% 3. PLOT 3D COMPLETO
% -------------------------------------------------------------------------
figure('Name','Trasferimento Orbitale - Soluzione Circolare');

ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal; view(35, 25);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Sequenza di Manovre e Nodi (Orbita Ausiliaria Circolare)','Color','w','FontSize',14);

% Terra con Mappa Topografica (Nativo in MATLAB)
load topo; 
[XE, YE, ZE] = sphere(50);
surf(XE*6371, YE*6371, ZE*6371, 'CData', topo, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'HandleVisibility', 'off');
colormap(ax, topomap1); % Applica i colori della terra

dth = deg2rad(1);

% --- PLOT ORBITE ---
% 1. Orbita Iniziale
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [0.2 0.6 1], '-', 1.5, '1. Orbita Iniziale');
% 2. Trasferimento Salita
r_a_t1 = a_aux; 
a_t1 = (r_p_i + r_a_t1)/2;
e_t1 = (r_a_t1 - r_p_i)/(r_a_t1 + r_p_i);
plotOrbit(a_t1, e_t1, i_i, OM_i, om_i, 0, pi, dth, mu, [1 0.6 0.1], '-', 2, '2. Trasferimento Salita');
% 3. Orbita Ausiliaria (pre-piano)
plotOrbit(a_aux, e_aux, i_i, OM_i, om_i, 0, 2*pi, dth, mu, [1 1 0], '--', 1.5, '3. Orbita Aux (Pre-Piano)');
% 4. Orbita Ausiliaria (post-piano)
plotOrbit(a_aux, e_aux, i_f, OM_f, om_tmp, 0, 2*pi, dth, mu, [0.8 0.2 0.8], '-.', 1.5, '4. Orbita Aux (Post-Piano)');
% 5. Trasferimento Discesa
a_t2 = (a_aux + r_p_f)/2;
e_t2 = (a_aux - r_p_f)/(a_aux + r_p_f);
plotOrbit(a_t2, e_t2, i_f, OM_f, om_f, pi, 2*pi, dth, mu, [1 0.3 0.3], ':', 2, '5. Trasferimento Discesa');
% 6. Orbita Finale Target
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, dth, mu, [0 1 0], '-', 2, '6. Orbita Target');

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

function r_ijk = getPos3D(a, e, i, OM, om, th)
    % GETPOS3D: Calcola il vettore posizione 3D [X; Y; Z] per una data anomalia vera
    r = (a*(1 - e^2)) / (1 + e*cos(th));
    
    % Vettore in componenti perifocali (piano dell'orbita)
    r_pqw = [r*cos(th); r*sin(th); 0];
    
    % Matrici di rotazione
    R3_OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
    R1_i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
    
    % Rotazione da perifocale a inerziale geocentrico (ECI)
    T_pqw2ijk = R3_OM * R1_i * R3_om;
    r_ijk = T_pqw2ijk * r_pqw;
end

%% =========================================================================
% 4. ANIMAZIONE 3D AVANZATA - TRACCIA MULTI-COLORE E ORBITE DI RIFERIMENTO
% =========================================================================
% Ricalcolo parametri ellissi di trasferimento (già calcolati per il plot 3D,
% li ridefiniamo per sicurezza e leggibilità della section)
r_p_i  = a_i*(1-e_i);
r_p_f  = a_f*(1-e_f);
a_t1   = (r_p_i + a_aux)/2;
e_t1   = (a_aux - r_p_i)/(a_aux + r_p_i);
a_t2   = (a_aux + r_p_f)/2;
e_t2   = (a_aux - r_p_f)/(a_aux + r_p_f);

% --- 1. SETUP SCENA ---
fig_anim = figure('Name', 'Simulazione Dinamica Trasferimento Orbitale Circolare', 'Color', 'w', 'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
hold on; grid on; axis equal; view(35, 25);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Animazione Sequenza Manovre (Orbita Ausiliaria Circolare)', 'FontSize', 14);

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
plotOrbitStaticAnim(a_i, e_i, i_i, OM_i, om_i, mu, [0.2 0.6 1], '--', 1, 'Rif. Iniziale');
plotOrbitStaticAnim(a_aux, e_aux, i_i, OM_i, om_i, mu, [1 1 0], '--', 1, 'Rif. Aux Pre-Piano');
plotOrbitStaticAnim(a_aux, e_aux, i_f, OM_f, om_f, mu, [0.8 0.2 0.8], '--', 1, 'Rif. Aux Post-Piano/Peri');
plotOrbitStaticAnim(a_f, e_f, i_f, OM_f, om_f, mu, [0 1 0], '--', 1, 'Rif. Finale Target');

% --- 3. DEFINIZIONE SEGMENTI E COLORI TRACCIA ---
% Ogni riga: {a, e, i, OM, om, th_start, th_end, nome, colore_traccia}
segmenti = {};

% Seg 1: Coasting su Orbita Iniziale (da th_i a pericentro th=0)
th_s = th_i; th_e = 0; if th_e <= th_s; th_e = th_e + 2*pi; end
segmenti{1} = {a_i, e_i, i_i, OM_i, om_i, th_s, th_e, 'Coasting Iniziale', [0.2 0.6 1]};

% Seg 2: Trasferimento in Salita (da pericentro iniziale th=0 ad apocentro th=pi)
segmenti{2} = {a_t1, e_t1, i_i, OM_i, om_i, 0, pi, 'Trasferimento Salita', [1 0.6 0.1]};

% Seg 3: Coasting su Aux Pre-Piano (da th=pi al nodo di cambio piano th_plane_tmp)
th_s = pi; th_e = th_plane_tmp; if th_e <= th_s; th_e = th_e + 2*pi; end
segmenti{3} = {a_aux, e_aux, i_i, OM_i, om_i, th_s, th_e, 'Verso Nodo Cambio Piano', [1 1 0]};

% --- FIX TELETRASPORTO (Usando la tua logica!) ---
% Il satellite si trova a th_plane_tmp. Ruotando il pericentro da om_tmp a om_f,
% la sua nuova anomalia vera è semplicemente scalata della differenza tra i due.
th_post_plane_new = mod(th_plane_tmp + om_tmp - om_f, 2*pi);

% Seg 4: Coasting su Aux Post-Piano/Peri (dal nodo al punto di discesa th=pi)
th_s = th_post_plane_new; 
th_e = pi; % Arrivo per la bitangente di discesa
if th_e <= th_s; th_e = th_e + 2*pi; end
segmenti{4} = {a_aux, e_aux, i_f, OM_f, om_f, th_s, th_e, 'Allineamento Discesa', [0.8 0.2 0.8]};

% Seg 5: Trasferimento in Discesa (da apocentro th=pi a pericentro th=2*pi)
segmenti{5} = {a_t2, e_t2, i_f, OM_f, om_f, pi, 2*pi, 'Trasferimento Discesa', [1 0.3 0.3]};

% Seg 6: Coasting su Orbita Finale Target (da pericentro th=0 a th_f)
th_s = 0; th_e = th_f; if th_e <= th_s; th_e = th_e + 2*pi; end
segmenti{6} = {a_f, e_f, i_f, OM_f, om_f, th_s, th_e, 'Arrivo a Target', [0 1 0]};

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
        r = getPosAnim(seg{1}, seg{2}, seg{3}, seg{4}, seg{5}, th_v(k), mu);
        
        % Aggiorna satellite
        set(h_sat, 'XData', r(1), 'YData', r(2), 'ZData', r(3));
        
        % Aggiorna traccia corrente
        t_x(end+1) = r(1); t_y(end+1) = r(2); t_z(end+1) = r(3);
        set(h_trail, 'XData', t_x, 'YData', t_y, 'ZData', t_z);
        
        drawnow;
        pause(0.005); % Regola velocità animazione
    end
    
    % Marker Manovra alla fine di ogni segmento
    plot3(t_x(end), t_y(end), t_z(end), 'x', 'MarkerEdgeColor', seg{9}, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
end

legend('show', 'Location', 'bestoutside', 'FontSize', 9);

% --- FUNZIONI INTERNE PER ANIMAZIONE ---
function r_ijk = getPosAnim(a, e, i, OM, om, th, mu)
    r_mag = (a*(1 - e^2)) / (1 + e*cos(th));
    r_pqw = [r_mag*cos(th); r_mag*sin(th); 0];
    R3_OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
    R1_i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
    r_ijk = R3_OM * R1_i * R3_om * r_pqw;
end

function plotOrbitStaticAnim(a, e, i, OM, om, mu, col, stile, width, nome)
    th = linspace(0, 2*pi, 200);
    pts = zeros(3, 200);
    for j = 1:200
        r_mag = (a*(1 - e^2)) / (1 + e*cos(th(j)));
        rpqw = [r_mag*cos(th(j)); r_mag*sin(th(j)); 0];
        R3OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
        R1i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
        R3om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
        pts(:,j) = R3OM * R1i * R3om * rpqw;
    end
    plot3(pts(1,:), pts(2,:), pts(3,:), 'Color', [col 0.4], 'LineStyle', stile, 'LineWidth', width, 'DisplayName', nome);
end