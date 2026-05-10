clear
close all
clc

% --- 1. CARICAMENTO DATI ---
T = load("DatiSC1-2026.txt");
mu = 398600.4;
gruppo = 23;

% Orbita finale (vettori di stato)
rr_f = T(gruppo,8:10)';
vv_f = T(gruppo,11:13)';
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

% Orbita iniziale (parametri kepleriani)
orbita_iv = T(gruppo,2:7);
a_i = orbita_iv(1); e_i = orbita_iv(2); i_i = orbita_iv(3);
OM_i = orbita_iv(4); om_i = orbita_iv(5); th_i = orbita_iv(6);


% --- 2. BIELLITTICA REALE (CON RITORNO LAMBERT) ---
ra_t = 170000; % [km] Apocentro lontano

% PUNTO 1: Partenza dal Pericentro Iniziale
rp_i = a_i * (1 - e_i);
[r1, v1_orb] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu); 

% Innalzamento verso ra_t
a_trans1 = (rp_i + ra_t) / 2;
e_trans1 = (ra_t - rp_i) / (ra_t + rp_i);
v1_trans_mag = sqrt(2*mu/rp_i - mu/a_trans1);
v1_trans1 = (v1_trans_mag / norm(v1_orb)) * v1_orb; 
dv1 = norm(v1_trans1 - v1_orb);

% PUNTO 2: Arrivo all'apocentro lontano (Punto reale nello spazio 3D)
[r2, v2_t1] = par2car(a_trans1, e_trans1, i_i, OM_i, om_i, pi, mu);

% PUNTO 3: Il nostro bersaglio esatto
r3 = rr_f;
v3_orb_f = vv_f;

% Usiamo Lambert per trovare l'orbita esatta di connessione tra r2 e r3
TOF_test = linspace(20000, 200000, 400); % Testiamo vari tempi di caduta
min_dv = inf;

for j = 1:length(TOF_test)
    tof = TOF_test(j);
    [v2_L, v3_L] = solve_lambert(r2, r3, tof, mu, 1);
    
    dv2_test = norm(v2_L - v2_t1);
    dv3_test = norm(v3_orb_f - v3_L);
    
    if (dv2_test + dv3_test) < min_dv
        min_dv = dv2_test + dv3_test;
        best_tof = tof;
        v2_lambert = v2_L;
        v3_lambert = v3_L;
        dv2 = dv2_test;
        dv3 = dv3_test;
    end
end

% --- 3. RISULTATI ---
DeltaV_Tot = dv1 + dv2 + dv3;
fprintf('\n--- RISULTATI BIELLITTICA REALE (Lambert Return) ---\n');
fprintf('DeltaV 1 (Innalzamento a Pericentro): %.4f km/s\n', dv1);
fprintf('DeltaV 2 (Correzione 3D a 150k km):   %.4f km/s\n', dv2);
fprintf('DeltaV 3 (Frenata al Pericentro Fin): %.4f km/s\n', dv3);
fprintf('DELTAV TOTALE:                        %.4f km/s\n\n', DeltaV_Tot);

% --- 4. TEMPI DI VOLO ---
dt_coast = TOF(a_i, e_i, th_i, 0, mu); 
dt_trans1 = pi * sqrt(a_trans1^3/mu);
dt_trans2 = best_tof; % Il tempo ottimale trovato da Lambert
TOF_Tot = dt_coast + dt_trans1 + dt_trans2;
fprintf('Tempo totale: %.2f giorni\n', TOF_Tot/(24*3600));

% --- 5. PLOT GRAFICO AVANZATO ---
figure('Name', 'Biellittica Geometricamente Esatta');
hold on; grid on; axis equal; view(3);

% A. Disegno della Terra
R_earth = 6371; 
[xE, yE, zE] = sphere(50); 
try
    load topo topo topomap1; 
    surf(xE * R_earth, yE * R_earth, zE * R_earth, 'FaceColor', 'texturemap', 'CData', topo, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    colormap(topomap1); 
catch
    surf(xE * R_earth, yE * R_earth, zE * R_earth, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
end

dth = 0.05;

% 1. Orbita Iniziale
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, dth, mu, '#0072BD', 'Orbita Iniziale');

% 2. Trasferimento 1 (Salita)
plotOrbit(a_trans1, e_trans1, i_i, OM_i, om_i, 0, pi, dth, mu, '#D95319', 'Trasferimento 1 (Salita)');

% 3. Trasferimento 2 (Discesa di Lambert)
% Ricaviamo i parametri della traiettoria di Lambert per poterla plottare
[a_L, e_L, i_L, OM_L, om_L, th2_L] = car2par(r2, v2_lambert, mu);
[~, ~, ~, ~, ~, th3_L] = car2par(r3, v3_lambert, mu);
if th3_L < th2_L
    th3_L = th3_L + 2*pi; % Gestione attraversamento asse
end
plotOrbit(a_L, e_L, i_L, OM_L, om_L, th2_L, th3_L, dth, mu, '#EDB120', 'Trasferimento 2 (Lambert)');

% 4. Orbita Finale
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, dth, mu, '#77AC30', 'Orbita Finale');

% C. Marker delle Manovre
plot3(r1(1), r1(2), r1(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Impulso 1');
plot3(r2(1), r2(2), r2(3), 's', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'DisplayName', 'Impulso 2');
plot3(r3(1), r3(2), r3(3), 'd', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'DisplayName', 'Impulso 3');

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title(sprintf('\\Delta V Totale Esatto = %.3f km/s', DeltaV_Tot));
legend('show', 'Location', 'best');
hold off;

% =========================================================================
function plotOrbit(a, e, i, OM, om, th0, thf, dth, mu, colore, nome_legenda)
    th_vec = th0:dth:thf;
    if th_vec(end) ~= thf, th_vec = [th_vec, thf]; end
    N = length(th_vec);
    rr = zeros(3, N);
    for k = 1:N
        [rr_i, ~] = par2car(a, e, i, OM, om, th_vec(k), mu);
        rr(:, k) = rr_i;
    end
    plot3(rr(1,:), rr(2,:), rr(3,:), 'Color', colore, 'LineWidth', 1.5, 'DisplayName', nome_legenda);
end