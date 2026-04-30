% =========================================================================
%   TRASFERIMENTO ORBITALE - SEQUENZA 5 MANOVRE
%   1) Partenza dal punto assegnato dell'orbita iniziale
%   2) Bitangente 'pa' verso orbita ausiliaria ampia  (ottimizzata su a_aux)
%   3) Cambio di piano sull'orbita ausiliaria (lontano dall'attrattore)
%   4) Bitangente per adattare le dimensioni all'orbita finale
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
gruppo = 23;

% --- Orbita target / finale (vettori di stato alla riga 23) --------------
rr_f   = T(gruppo, 8:10)';
vv_f   = T(gruppo, 11:13)';

[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

fprintf('=== ORBITA FINALE ===\n');
fprintf('  a  = %.4f km\n', a_f);
fprintf('  e  = %.6f\n',    e_f);
fprintf('  i  = %.4f deg\n', rad2deg(i_f));
fprintf('  OM = %.4f deg\n', rad2deg(OM_f));
fprintf('  om = %.4f deg\n', rad2deg(om_f));
fprintf('  th = %.4f deg\n', rad2deg(th_f));

% -------------------------------------------------------------------------
% 2. DEFINIZIONE ORBITA INIZIALE
%    (modifica questi parametri secondo il tuo caso specifico)
% -------------------------------------------------------------------------

orbita_iv = T(gruppo,2:7);
a_i  = orbita_iv(1);
e_i  = orbita_iv(2);
i_i  = orbita_iv(3);
OM_i = orbita_iv(4);
om_i = orbita_iv(5);
th_i = orbita_iv(6);

fprintf('\n=== ORBITA INIZIALE ===\n');
fprintf('  a  = %.4f km\n', a_i);
fprintf('  e  = %.6f\n',    e_i);
fprintf('  i  = %.4f deg\n', rad2deg(i_i));
fprintf('  OM = %.4f deg\n', rad2deg(OM_i));
fprintf('  om = %.4f deg\n', rad2deg(om_i));
fprintf('  th = %.4f deg\n', rad2deg(th_i));

% -------------------------------------------------------------------------
% 3. OTTIMIZZAZIONE ORBITA AUSILIARIA (ciclo for su a_aux)
%    La bitangente 'pa' parte dal pericentro dell'orbita iniziale e
%    arriva all'apocentro dell'orbita ausiliaria.
%    Vogliamo massimizzare a_aux riducendo al minimo il costo TOTALE
%    (manovra 2 + 3 + 4), oppure semplicemente trovare il minimo del
%    costo complessivo della sequenza rispetto ad a_aux.
%
%    L'orbita ausiliaria è costruita come segue:
%      - r_p_aux = r_p_i  (stessa quota al pericentro della iniziale)
%      - a_aux   variabile → e_aux = 1 - r_p_aux / a_aux
%    Dopo la bitangente 'pa' si è all'apocentro di aux.
%    Il cambio di piano avviene lì (punto più lontano → velocità minima).
%    Poi si esegue una bitangente 'ap' per tornare alle dimensioni
%    dell'orbita finale.
% -------------------------------------------------------------------------

r_p_i = a_i*(1-e_i);                    % raggio pericentro orbita iniziale
r_a_f = a_f*(1+e_f);                    % raggio apocentro orbita finale

% Range di a_aux da esplorare
%   - minimo: a tale che r_p_aux = r_p_i e r_a_aux >= r_a_f
%   - massimo: orbita molto ampia (es. 10x la finale)

a_aux_min = (r_p_i + r_a_f) / 2;        % la più piccola possibile (bitangente diretta pa→finale)
a_aux_max = 800000;
N_search  = 5000;
a_aux_vec = linspace(a_aux_min, a_aux_max, N_search);

dv_total_vec = zeros(1, N_search);

for k = 1:N_search
    a_aux_k = a_aux_vec(k);
    e_aux_k = 1 - r_p_i / a_aux_k;     % r_p_aux = r_p_i  → bitangente 'pa'

    % Se l'eccentricità supera 1 (iperbolica) saltiamo
    if e_aux_k >= 1 || e_aux_k < 0
        dv_total_vec(k) = Inf;
        continue;
    end

    % -- Manovra 2: bitangente 'pa' (pericentro i → apocentro aux) --------
    [dv2_1, dv2_2, ~] = bitangentTransfer(a_i, e_i, a_aux_k, e_aux_k, 'pa', mu);

    % -- Manovra 3: cambio di piano all'apocentro dell'orbita aux ---------
    %    Usiamo la funzione changeOrbitalPlane sull'orbita aux.
    %    Dopo la manovra 2 siamo all'apocentro aux, quindi:
    %      - i_i, OM_i dell'orbita aux = i dell'orbita iniziale, OM_i
    %      - vogliamo i_f, OM_f dell'orbita finale
    try
        [dv3, ~, ~] = changeOrbitalPlane(a_aux_k, e_aux_k, i_i, OM_i, om_i, i_f, OM_f, mu);
    catch
        dv_total_vec(k) = Inf;
        continue;
    end

    % -- Manovra 4: bitangente per passare da aux a orbita con stesse
    %    dimensioni della finale.
    %    Siamo ancora all'apocentro dell'orbita aux.
    %    Vogliamo arrivare a un'orbita con a = a_f, e = e_f.
    %    Usiamo tipo 'ap': apocentro aux → pericentro finale
    r_a_aux = a_aux_k*(1+e_aux_k);
    r_p_f   = a_f*(1-e_f);
    if r_a_aux < r_p_f
        dv_total_vec(k) = Inf;
        continue;
    end
    [dv4_1, dv4_2, ~] = bitangentTransfer(a_aux_k, e_aux_k, a_f, e_f, 'ap', mu);

    dv_total_vec(k) = abs(dv2_1) + abs(dv2_2) + abs(dv3) + abs(dv4_1) + abs(dv4_2);
end

% Trova il minimo
[~, idx_opt] = min(dv_total_vec);
a_aux = a_aux_vec(idx_opt);
e_aux = 1 - r_p_i / a_aux;

fprintf('\n=== OTTIMIZZAZIONE ORBITA AUSILIARIA ===\n');
fprintf('  a_aux ottimale = %.4f km\n', a_aux);
fprintf('  e_aux          = %.6f\n',    e_aux);
fprintf('  DeltaV totale (manovre 2-4) = %.6f km/s\n', dv_total_vec(idx_opt));

% -------------------------------------------------------------------------
% 4. CALCOLO DETTAGLIATO DI TUTTE LE MANOVRE
% -------------------------------------------------------------------------
% Parametri aggiornati lungo la sequenza
% -- Orbita corrente: iniziale
a_cur  = a_i;
e_cur  = e_i;
i_cur  = i_i;
OM_cur = OM_i;
om_cur = om_i;
th_cur = th_i;

fprintf('\n========================================================\n');
fprintf('   DETTAGLIO MANOVRE\n');
fprintf('========================================================\n');

% ===== COASTING 1: dal punto di partenza al pericentro ==================
%  La bitangente 'pa' parte dal pericentro → coasting da th_i a th=0
th_peri = 0;   % anomalia vera al pericentro
dt_coast1 = TOF(a_cur, e_cur, th_cur, th_peri, mu);
if dt_coast1 < 0   % se già oltre il pericentro, giro completo
    T_orb = 2*pi*sqrt(a_cur^3/mu);
    dt_coast1 = dt_coast1 + T_orb;
end
fprintf('\n--- COASTING 1: th=%.1f deg → pericentro (th=0) ---\n', rad2deg(th_cur));
fprintf('    Durata = %.2f s  (%.4f ore)\n', dt_coast1, dt_coast1/3600);

% ===== MANOVRA 2: Bitangente 'pa' (pericentro i → apocentro aux) ========
[dv2_1, dv2_2, dt_m2] = bitangentTransfer(a_i, e_i, a_aux, e_aux, 'pa', mu);
fprintf('\n--- MANOVRA 2: Bitangente ''pa'' i → aux ---\n');
fprintf('    DeltaV1 = %.6f km/s  (al pericentro di i)\n', dv2_1);
fprintf('    DeltaV2 = %.6f km/s  (all''apocentro di aux)\n', dv2_2);
fprintf('    TOF     = %.2f s  (%.4f ore)\n', dt_m2, dt_m2/3600);
fprintf('    |DeltaV_tot_m2| = %.6f km/s\n', abs(dv2_1)+abs(dv2_2));

% Aggiorno i parametri: ora sono sull'orbita aux, all'apocentro
a_cur  = a_aux;
e_cur  = e_aux;
% i, OM, om rimangono uguali all'orbita iniziale dopo bitangente
th_cur = pi;    % apocentro

% ===== MANOVRA 3: Cambio di piano (all'apocentro di aux) ================
[dv3, om_after_plane, th_plane] = changeOrbitalPlane(a_cur, e_cur, ...
                                   i_i, OM_i, om_i, ...
                                   i_f, OM_f, mu);

fprintf('\n--- MANOVRA 3: Cambio di piano ---\n');
fprintf('    DeltaV  = %.6f km/s\n', dv3);
fprintf('    Anomalia vera alla manovra = %.4f deg\n', rad2deg(th_plane));
fprintf('    om dopo cambio piano       = %.4f deg\n', rad2deg(om_after_plane));
fprintf('    i: %.4f → %.4f deg\n', rad2deg(i_i), rad2deg(i_f));
fprintf('    OM: %.4f → %.4f deg\n', rad2deg(OM_i), rad2deg(OM_f));

% Aggiorno parametri
i_cur  = i_f;
OM_cur = OM_f;
om_cur = om_after_plane;
th_cur = th_plane;

% ===== COASTING 2: dal punto del cambio piano all'apocentro =============
%  La bitangente 'ap' parte dall'apocentro di aux (th=pi)
%  Verifico se siamo già all'apocentro o dobbiamo attendere
th_apo = pi;
dt_coast2 = TOF(a_cur, e_cur, th_cur, th_apo, mu);
if dt_coast2 < 0
    T_orb = 2*pi*sqrt(a_cur^3/mu);
    dt_coast2 = dt_coast2 + T_orb;
end
fprintf('\n--- COASTING 2: th=%.1f deg → apocentro (th=180) ---\n', rad2deg(th_cur));
fprintf('    Durata = %.2f s  (%.4f ore)\n', dt_coast2, dt_coast2/3600);
th_cur = th_apo;

% ===== MANOVRA 4: Bitangente 'ap' (apocentro aux → finale) ==============
[dv4_1, dv4_2, dt_m4] = bitangentTransfer(a_aux, e_aux, a_f, e_f, 'ap', mu);
fprintf('\n--- MANOVRA 4: Bitangente ''ap'' aux → finale ---\n');
fprintf('    DeltaV1 = %.6f km/s  (all''apocentro di aux)\n', dv4_1);
fprintf('    DeltaV2 = %.6f km/s  (al pericentro di f)\n', dv4_2);
fprintf('    TOF     = %.2f s  (%.4f ore)\n', dt_m4, dt_m4/3600);
fprintf('    |DeltaV_tot_m4| = %.6f km/s\n', abs(dv4_1)+abs(dv4_2));

% Aggiorno: sono sul pericentro dell'orbita finale (a_f, e_f, i_f, OM_f)
% ma con om_cur (che può differire da om_f)
a_cur  = a_f;
e_cur  = e_f;
th_cur = 0;    % pericentro

% ===== MANOVRA 5: Rotazione argomento del pericentro ====================
[dv5, thi_vec, thf_vec] = changePericenterArg(a_cur, e_cur, om_cur, om_f, mu);

% Scegliamo il punto di manovra più conveniente (minima velocità → apocentro)
% thi ha due soluzioni: scegliamo quella all'apocentro (th vicino a pi)
[~, idx5] = min(abs(thi_vec - pi));
th_manov5 = thi_vec(idx5);
th_after5 = thf_vec(idx5);

fprintf('\n--- MANOVRA 5: Rotazione argomento pericentro ---\n');
fprintf('    om: %.4f → %.4f deg\n', rad2deg(om_cur), rad2deg(om_f));
fprintf('    DeltaV = %.6f km/s\n', dv5);
fprintf('    Anomalia vera alla manovra = %.4f deg\n', rad2deg(th_manov5));
fprintf('    Anomalia vera dopo manovra = %.4f deg\n', rad2deg(th_after5));

% ===== COASTING 3: dal pericentro al punto di manovra 5 =================
dt_coast3 = TOF(a_cur, e_cur, th_cur, th_manov5, mu);
if dt_coast3 < 0
    T_orb = 2*pi*sqrt(a_cur^3/mu);
    dt_coast3 = dt_coast3 + T_orb;
end
fprintf('\n--- COASTING 3: th=%.1f deg → th=%.1f deg (prima di manovra 5) ---\n', ...
        rad2deg(th_cur), rad2deg(th_manov5));
fprintf('    Durata = %.2f s  (%.4f ore)\n', dt_coast3, dt_coast3/3600);

% Aggiorno parametri finali
om_cur = om_f;
th_cur = th_after5;

% ===== RIEPILOGO BUDGET DeltaV ==========================================
dv_m2  = abs(dv2_1) + abs(dv2_2);
dv_m3  = abs(dv3);
dv_m4  = abs(dv4_1) + abs(dv4_2);
dv_m5  = abs(dv5);
dv_TOT = dv_m2 + dv_m3 + dv_m4 + dv_m5;

fprintf('\n========================================================\n');
fprintf('   BUDGET DELTA-V\n');
fprintf('========================================================\n');
fprintf('  Manovra 2 (bitangente pa)      : %.6f km/s\n', dv_m2);
fprintf('  Manovra 3 (cambio piano)        : %.6f km/s\n', dv_m3);
fprintf('  Manovra 4 (bitangente ap)       : %.6f km/s\n', dv_m4);
fprintf('  Manovra 5 (rot. pericentro)     : %.6f km/s\n', dv_m5);
fprintf('  -----------------------------------------------\n');
fprintf('  TOTALE                          : %.6f km/s\n', dv_TOT);
fprintf('========================================================\n');

dt_TOT = dt_coast1 + dt_m2 + dt_coast2 + dt_m4 + dt_coast3;
fprintf('\n  Tempo totale (coasting+trasferimenti): %.2f s  (%.4f h  = %.4f giorni)\n', ...
        dt_TOT, dt_TOT/3600, dt_TOT/86400);

% =========================================================================
% 5.  PLOT 3D
% =========================================================================
figure('Name','Trasferimento Orbitale - Vista 3D','Color','k','Position',[100 100 1200 900]);
ax = axes;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w');
hold on; grid on; axis equal;
xlabel('X [km]','Color','w'); ylabel('Y [km]','Color','w'); zlabel('Z [km]','Color','w');
title('Sequenza di Manovre Orbitali','Color','w','FontSize',14);

dth = deg2rad(0.5);   % passo angolare per il plot

% -- Orbita iniziale (intera) ------------------------------------------------
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi-dth, dth, mu);
h_init = findobj(gca,'Type','Line');
set(h_init(1),'Color',[0.2 0.6 1],'DisplayName','Orbita Iniziale');

% -- Coasting 1: da th_i a pericentro ----------------------------------------
if th_i ~= 0
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, th_i, 0, dth, mu);
    h_c1 = findobj(gca,'Type','Line');
    set(h_c1(1),'Color',[0.5 0.5 0.5],'LineStyle','--','DisplayName','Coasting 1');
end

% -- Trasferimento bitangente 1 (orbita di trasferimento pa) ------------------
% L'orbita di trasferimento ha: a_t2, e_t2, i_i, OM_i, om_i  (pericentro=r_pi, apocentro=r_a_aux)
r_p_i2 = a_i*(1-e_i);
r_a_aux2= a_aux*(1+e_aux);
a_t2   = (r_p_i2 + r_a_aux2)/2;
e_t2   = (r_a_aux2 - r_p_i2)/(r_a_aux2 + r_p_i2);
plotOrbit(a_t2, e_t2, i_i, OM_i, om_i, 0, pi, dth, mu);
h_t2 = findobj(gca,'Type','Line');
set(h_t2(1),'Color',[1 0.6 0.1],'LineStyle','-','DisplayName','Trasferimento 2 (bitangente pa)');

% -- Orbita ausiliaria (intera) -----------------------------------------------
plotOrbit(a_aux, e_aux, i_i, OM_i, om_i, 0, 2*pi-dth, dth, mu);
h_aux1 = findobj(gca,'Type','Line');
set(h_aux1(1),'Color',[1 1 0],'LineStyle','-','DisplayName','Orbita Aux (pre-piano)');

% -- Orbita ausiliaria dopo cambio di piano -----------------------------------
plotOrbit(a_aux, e_aux, i_f, OM_f, om_after_plane, 0, 2*pi-dth, dth, mu);
h_aux2 = findobj(gca,'Type','Line');
set(h_aux2(1),'Color',[0.8 0.2 0.8],'LineStyle','-','DisplayName','Orbita Aux (post-piano)');

% -- Trasferimento bitangente 2 (ap: apocentro aux → pericentro finale) -------
r_a_aux3 = a_aux*(1+e_aux);
r_p_f2   = a_f*(1-e_f);
a_t4     = (r_a_aux3 + r_p_f2)/2;
e_t4     = (r_a_aux3 - r_p_f2)/(r_a_aux3 + r_p_f2);
plotOrbit(a_t4, e_t4, i_f, OM_f, om_after_plane, pi, 2*pi, dth, mu);
h_t4 = findobj(gca,'Type','Line');
set(h_t4(1),'Color',[1 0.3 0.3],'LineStyle','-','DisplayName','Trasferimento 4 (bitangente ap)');

% -- Orbita intermedia (a_f, e_f, i_f, OM_f, om_after_plane) PRIMA di rot. om -
plotOrbit(a_f, e_f, i_f, OM_f, om_after_plane, 0, 2*pi-dth, dth, mu);
h_int = findobj(gca,'Type','Line');
set(h_int(1),'Color',[0.2 0.9 0.5],'LineStyle','--','DisplayName','Orbita intermedia (pre-rot)');

% -- Orbita finale (intera) --------------------------------------------------
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi-dth, dth, mu);
h_fin = findobj(gca,'Type','Line');
set(h_fin(1),'Color',[0.1 1 0.1],'DisplayName','Orbita Finale');

% -- Punto di partenza -------------------------------------------------------
[rr_start, ~] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i, mu);
scatter3(rr_start(1), rr_start(2), rr_start(3), 80, 'c', 'filled', ...
         'DisplayName','Punto di Partenza');

% -- Punto di manovra 3 (cambio piano) --------------------------------------
[rr_mp3, ~] = par2car(a_aux, e_aux, i_i, OM_i, om_i, th_plane, mu);
scatter3(rr_mp3(1), rr_mp3(2), rr_mp3(3), 80, [1 0.6 0], 'filled', ...
         'DisplayName','Manovra 3 (cambio piano)');

% -- Punto di manovra 5 (rot. pericentro) ------------------------------------
[rr_mp5, ~] = par2car(a_f, e_f, i_f, OM_f, om_after_plane, th_manov5, mu);
scatter3(rr_mp5(1), rr_mp5(2), rr_mp5(3), 80, [1 0.3 0.3], 'filled', ...
         'DisplayName','Manovra 5 (rot. pericentro)');

% -- Punto di arrivo (orbita finale) -----------------------------------------
[rr_end, ~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);
scatter3(rr_end(1), rr_end(2), rr_end(3), 100, [0 1 0], 'p', 'filled', ...
         'DisplayName','Punto di Arrivo');

legend('show','Location','bestoutside','TextColor','w','Color','k','FontSize',8);
view(35, 25);

% =========================================================================
% 6.  PLOT CURVA DI OTTIMIZZAZIONE
% =========================================================================
figure('Name','Ottimizzazione a_aux','Color','w');
plot(a_aux_vec/1e3, dv_total_vec, 'b-', 'LineWidth', 1.5); hold on;
xline(a_aux/1e3, 'r--', 'LineWidth', 1.5);
yline(dv_total_vec(idx_opt), 'g--', 'LineWidth', 1.2);
scatter(a_aux/1e3, dv_total_vec(idx_opt), 100, 'r', 'filled');
xlabel('a_{aux} [×10^3 km]'); ylabel('\DeltaV totale [km/s]');
title('Ottimizzazione semi-asse orbita ausiliaria');
grid on; legend({'DeltaV totale', 'a_{aux} ottimale'}, 'Location','best');

fprintf('\nScript completato con successo.\n');