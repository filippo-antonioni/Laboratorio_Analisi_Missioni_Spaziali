clear
clc
close all

%% SECTION 1 - Geocentrica -> Eliocentrica (Iperbole di Fuga)
%%% Load dati orbita finale geocentrica (scenario 1) - ORBITA DI PARCHEGGIO %%%
T = load("DatiSC1-2026.txt");
mu_earth = 398600;
m_earth = 5.974e24; % [kg]
gruppo = 23;

% Orbita target/finale (vettori di stato)
rr_f = T(gruppo,8:10)';
vv_f = T(gruppo,11:13)';
state_f = [rr_f; vv_f];

% Trasformazione per avere i parametri kepleriani bersaglio (corretto: mu_earth)
[a_op, e_op, i_op, OM_op, om_op, th_op] = car2par(rr_f, vv_f, mu_earth);

% Moduli del raggio e velocità di pericentro dell'orbita di parcheggio
r_op = a_op * (1 - e_op);
v_op = sqrt(mu_earth * (2/r_op - 1/a_op));

%%% Load dati orbita eliocentrica (scenario 2) %%%
load('workspace_sc2_gridsearch.mat'); 
mu_sun  = 1.32712440018e11; % Parametro gravitazionale del Sole [km^3/s^2]
m_sun   = 1.989e30;         % [kg]
R_sun   = 696340;           % Raggio del Sole [km]
AU      = 149597870.7;      % 1 Unità Astronomica [km]

%%% Calcolo SOI terra (usando la distanza reale alla partenza) %%%
% r1_opt è stato calcolato nello Scenario 2 (o lo ricalcoli con ottimo_th1)
R_terra_partenza = norm(r1_opt); 
r_SOI_terra = R_terra_partenza * (m_earth/m_sun)^(2/5);

% Velocità eliocentrica della Terra alla partenza (V1i)
[r1i, V1i] = par2car(a_T, e_T, i_T, OM_T, om_T, ottimo_th1, mu_sun);

% Velocità eliocentrica del satellite sull'orbita di trasferimento alla partenza (V1T)
[R1T, V1T] = par2car(ottimo_aT, ottimo_eT, i_trasf_opt, OM_transf_opt, ottimo_omT, th1_T_opt, mu_sun);

% Eccesso iperbolico vettoriale e modulo
v_inf_1_vec = V1T - V1i;
v_inf_1 = norm(v_inf_1_vec);

%%% 2. CARATTERIZZAZIONE IPERBOLE DI FUGA %%%
r_PH1 = r_op;                    % Impongo il raggio di pericentro 
a_H1 = - mu_earth / (v_inf_1^2); % Semiasse maggiore dell'iperbole (< 0)
e_H1 = 1 - (r_PH1 / a_H1);       % Eccentricità dell'iperbole (> 1)

%%% 3. CALCOLO DELTA V DI FUGA %%%
% Velocità al pericentro dell'iperbole
v_PH1 = sqrt(v_inf_1^2 + (2*mu_earth)/r_PH1); 
% Delta V richiesto per la manovra
DeltaV_1 = abs(v_PH1 - v_op);

%%% 4. CONTROLLO DI COERENZA DELL'IPERBOLE %%%
if e_H1 > 1 && a_H1 < 0
    fprintf('\n[OK] Controllo superato: L''orbita di fuga è un''iperbole (e = %.4f, a = %.2f km).\n', e_H1, a_H1);
else
    warning('ATTENZIONE: I parametri calcolati non corrispondono a un''iperbole!');
end

%%% 5. PLOT DEL SISTEMA GEOCENTRICO 3D (Riscalato e con Terra realistica) %%%
figure('Name', 'Scenario 3: Iperbole di Fuga dalla Terra (3D)');
hold on; grid on; axis equal; view(3);

% Disegno la Terra con Topomap realistica integrata in MATLAB
R_earth = 6371; 
[X_E, Y_E, Z_E] = sphere(50);
load topo; % Carica dataset topografico nativo
surf(X_E*R_earth, Y_E*R_earth, Z_E*R_earth, 'CData', topo, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', 'DisplayName', 'Terra');
colormap(topomap1); % Applica la mappa colori terrestre

% Definizione vettori anomalie per i plot
th_ell = linspace(0, 2*pi, 200);
th_inf = acos(-1/e_H1); 
th_hyp = linspace(-th_inf + 0.05, th_inf - 0.05, 200);

% A. Orbita di parcheggio (Ellisse - Scenario 1)
r_op_plot = zeros(3, length(th_ell));
for k = 1:length(th_ell)
    [r_op_plot(:,k), ~] = par2car(a_op, e_op, i_op, OM_op, om_op, th_ell(k), mu_earth);
end
plot3(r_op_plot(1,:), r_op_plot(2,:), r_op_plot(3,:), 'b', 'LineWidth', 1.5, 'DisplayName', 'Orbita di Parcheggio');

% B. Iperbole di fuga
r_hyp_plot = zeros(3, length(th_hyp));
for k = 1:length(th_hyp)
    [r_hyp_plot(:,k), ~] = par2car(a_H1, e_H1, i_op, OM_op, om_op, th_hyp(k), mu_earth);
end
plot3(r_hyp_plot(1,:), r_hyp_plot(2,:), r_hyp_plot(3,:), 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Fuga');

% C. Marker nel punto di manovra (Pericentro)
[r_manovra, ~] = par2car(a_op, e_op, i_op, OM_op, om_op, 0, mu_earth);
plot3(r_manovra(1), r_manovra(2), r_manovra(3), 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV: Punto di Manovra');

% Riscalamento assi in base alle dimensioni dell'orbita di parcheggio
r_ap_op = a_op * (1 + e_op); % Raggio di apocentro
lim_zoom = r_ap_op * 2;      % Finestra ampia il doppio dell'apocentro
xlim([-lim_zoom lim_zoom]); ylim([-lim_zoom lim_zoom]); zlim([-lim_zoom lim_zoom]);

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Iniezione su Iperbole di Fuga Geocentrica (Vista 3D)');
legend('Location', 'best');
hold off;


%%% 6. PLOT DEL SISTEMA GEOCENTRICO 2D (Piano Orbitale) %%%
% Nel piano perifocale z = 0, usiamo equazioni polari per disegnare le curve
figure('Name', 'Scenario 3: Iperbole di Fuga (2D Piano Orbitale)');
hold on; grid on; axis equal;

% Terra in 2D con Topomap (Generiamo la sfera ma la guarderemo dall'alto)
load topo;
[X_E2, Y_E2, Z_E2] = sphere(50);
surf(X_E2*R_earth, Y_E2*R_earth, Z_E2*R_earth, 'CData', topo, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', 'DisplayName', 'Terra');
colormap(topomap1);

% A. Orbita di Parcheggio 2D
r_op_mag = a_op * (1 - e_op^2) ./ (1 + e_op * cos(th_ell));
% Creiamo un array Z pari al raggio terrestre per far "galleggiare" la linea sopra la sfera
z_layer_ell = ones(size(th_ell)) * R_earth; 
plot3(r_op_mag .* cos(th_ell), r_op_mag .* sin(th_ell), z_layer_ell, 'b', 'LineWidth', 1.5, 'DisplayName', 'Orbita di Parcheggio');

% B. Iperbole di Fuga 2D
r_hyp_mag = a_H1 * (1 - e_H1^2) ./ (1 + e_H1 * cos(th_hyp));
z_layer_hyp = ones(size(th_hyp)) * R_earth;
plot3(r_hyp_mag .* cos(th_hyp), r_hyp_mag .* sin(th_hyp), z_layer_hyp, 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Fuga');

% C. Marker Manovra 2D (Il pericentro è per definizione a th = 0, sull'asse X positivo)
plot3(r_op, 0, R_earth, 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV: Punto di Manovra');

% Forziamo la vista 2D (perfettamente dall'alto)
view(2);

% Formattazione grafico 2D
xlim([-lim_zoom lim_zoom]); ylim([-lim_zoom lim_zoom]);
xlabel('Asse Pericentrale P [km]'); ylabel('Asse Trasverso Q [km]');
title('Vista 2D sul Piano Orbitale (Sistema Perifocale)');
legend('Location', 'best');
hold off;


%% SECTION 2 - Eliocentrica -> Asteroide

% Load dati orbita dell'asteroide 363505 (2003 UC20)
a_A  = 0.781241 * AU;       % Semiasse maggiore [km]
e_A  = 0.336932;            % Eccentricità [-]
i_A  = 3.78 * (pi/180);     % Inclinazione [rad]
OM_A = 187.92 * (pi/180);   % Anomalia del nodo ascendente (RAAN) [rad]
om_A = 60.16 * (pi/180);    % Argomento del pericentro [rad]
m_ast = 6.914e12; % [kg]
diam_ast = 1.88; % [km]

%%% Calcolo SOI asteroide (usando la distanza reale all'arrivo) %%%
R_ast_arrivo = norm(r2_opt);
r_SOI_ast = R_ast_arrivo * (m_ast/m_sun)^(2/5);

% Velocità eliocentrica dell'asteroide all'arrivo (V2f)
[r2f, V2f] = par2car(a_A, e_A, i_A, OM_A, om_A, ottimo_th2, mu_sun);

% Velocità eliocentrica del satellite sull'orbita di trasferimento all'arrivo (V2T)
[R2T, V2T] = par2car(ottimo_aT, ottimo_eT, i_trasf_opt, OM_transf_opt, ottimo_omT, th2_T_opt, mu_sun);

% Eccesso iperbolico vettoriale e modulo
v_inf_2_vec = V2f - V2T;
v_inf_2 = norm(v_inf_2_vec);

%%% CARATTERIZZAZIONE IPERBOLE DI ARRIVO %%%
G = 6.67430e-11; % [m^3 / (kg * s^2)] % Costante di gravitazione universale
mu_ast = G * m_ast / (1e9);
a_H2 = - mu_ast / v_inf_2^2;

%%% CONTROLLO SEMIASSE MAGGIORE DELL'IPERBOLE %%%
if a_H2 < 0
    fprintf('\n[OK] Controllo superato: L''orbita di fuga è un''iperbole (a = %.2f km).\n', a_H2);
else
    error('ATTENZIONE: I parametri calcolati non corrispondono a un''iperbole!');
end

% Creazione vettori di scansione
rvec = linspace(diam_ast/2, r_SOI_ast, 100);
e_vec = linspace(0, 1, 1000);
e_vec = e_vec(2:end-1); % Rimuovo 0 (è il caso circolare) e 1 (parabola)

N_r = length(rvec);
N_e = length(e_vec);

% --- PRE-ALLOCAZIONE MEMORIA (Per velocità e ordine) ---
% Struct per orbite circolari (vettori 1D)
orb_circ.r_ph  = zeros(1, N_r);
orb_circ.v     = zeros(1, N_r);
orb_circ.e_h   = zeros(1, N_r);
orb_circ.delta = zeros(1, N_r);
orb_circ.dv    = zeros(1, N_r);

% Struct per orbite ellittiche (matrici 2D: righe = rvec, colonne = e_vec)
orb_ell.r_a   = zeros(N_r, N_e);
orb_ell.e     = zeros(N_r, N_e);
orb_ell.r_p   = zeros(N_r, N_e);
orb_ell.a     = zeros(N_r, N_e);
orb_ell.e_h   = zeros(N_r, N_e);
orb_ell.delta = zeros(N_r, N_e);
orb_ell.dv    = zeros(N_r, N_e);

% Disabilitiamo temporaneamente i warning a schermo dentro il loop
% per evitare che la console esploda con 100.000 messaggi
fprintf('\nCalcolo sweep parametri in corso...\n');

for i = 1:N_r
    r_k = rvec(i);
    
    % ==========================================
    % 1. CASO ORBITA CIRCOLARE
    % ==========================================
    r_ph_k = r_k; 
    v_k_circ = sqrt(mu_ast / r_ph_k);
    e_h_k = 1 - (r_ph_k / a_H2);
    v_fuga_k_circ = sqrt(2 * mu_ast / r_ph_k);
    v_ph_k = sqrt(v_inf_2^2 + v_fuga_k_circ^2);
    
    % Salvataggio parametri circolari
    orb_circ.r_ph(i) = r_ph_k;
    orb_circ.v(i)    = v_k_circ;
    orb_circ.e_h(i)  = e_h_k;
    
    if e_h_k > 1
        orb_circ.delta(i) = -a_H2 * sqrt(e_h_k^2 - 1);
        orb_circ.dv(i)    = abs(v_k_circ - v_ph_k);
    else
        % Se non è un'iperbole, assegno infinito/NaN
        orb_circ.delta(i) = NaN;
        orb_circ.dv(i)    = Inf;        
    end
    
    % ==========================================
    % 2. CASO ORBITA ELLITTICA (Ciclo Annidato)
    % ==========================================
    
    %%% Caso orbita ellittica %%%
    r_a_k = r_k;
    
    for j = 1:N_e
        e_k = e_vec(j);
        
        % Calcolo subito il pericentro
        r_p_k = r_a_k * (1 - e_k) / (1 + e_k);
        
        % --- CHIAMO LA NOSTRA FUNZIONE DI CONTROLLO ---
        % Passiamo il pericentro, l'apocentro, il raggio dell'asteroide, la SOI e 100m di margine
        [orbita_ok, motivo] = check_feasibility(r_p_k, r_a_k, diam_ast/2, 0.1);
        
        if orbita_ok
            % Se l'orbita è geometricamente fattibile, calcoliamo la fisica
            r_ph_ell = r_p_k; 
            e_h_ell = 1 - (r_ph_ell / a_H2);
            
            if e_h_ell > 1
                v_fuga_k = sqrt(2 * mu_ast / r_ph_ell);
                v_ph_ell = sqrt(v_inf_2^2 + v_fuga_k^2);
                
                a_k = r_a_k / (1 + e_k);
                v_p_k = sqrt(mu_ast * (2/r_ph_ell - 1/a_k));
                
                % Salvataggio Dati Reali
                orb_ell.r_a(i, j)   = r_a_k;
                orb_ell.e(i, j)     = e_k;
                orb_ell.r_p(i, j)   = r_p_k;
                orb_ell.dv(i, j)    = abs(v_p_k - v_ph_ell);
                orb_ell.delta(i, j) = -a_H2 * sqrt(e_h_ell^2 - 1);
            else
                % Scarto: Non è un'iperbole
                orb_ell.dv(i, j) = NaN;
            end
        else
            % Scarto: L'orbita NON ha passato il controllo geometrico (Impatto, Fuga, ecc.)
            orb_ell.dv(i, j) = NaN;
            
            % Se volessi fare debug per capire quante schiantano, potresti usare disp(motivo), 
            % ma sconsiglio di scommentarlo altrimenti ti intasa la console.
            % disp(motivo);
        end
    end
end
fprintf('Sweep completato con successo!\n');

% =========================================================================
% --- ESTRAZIONE DEI RISULTATI OTTIMALI ---
% =========================================================================

% 1. Ottimo per Orbita Circolare
[dv_circ_min, idx_circ] = min(orb_circ.dv);
r_circ_opt = orb_circ.r_ph(idx_circ);
e_h_circ_opt = orb_circ.e_h(idx_circ);

% 2. Ottimo per Orbita Ellittica
% min(matrice(:)) trova il minimo assoluto in tutta la matrice 2D
[dv_ell_min, idx_ell_linear] = min(orb_ell.dv(:));
% ind2sub converte l'indice lineare in riga (raggio) e colonna (eccentricità)
[row_ell, col_ell] = ind2sub(size(orb_ell.dv), idx_ell_linear);

a_ell_opt   = orb_ell.a(row_ell, col_ell);
e_ell_opt   = orb_ell.e(row_ell, col_ell);
r_p_ell_opt = orb_ell.r_p(row_ell, col_ell);
r_a_ell_opt = orb_ell.r_a(row_ell, col_ell);
e_h_ell_opt = orb_ell.e_h(row_ell, col_ell);

% =========================================================================
% --- STAMPA DEI RISULTATI A SCHERMO ---
% =========================================================================

fprintf('\n======================================================\n');
fprintf('--- ANALISI MANOVRA DI CATTURA (2003 UC20) ---\n');
fprintf('======================================================\n');

fprintf('\n>>> 1. MIGLIOR ORBITA CIRCOLARE\n');
fprintf('Raggio Orbita  : %.4f km\n', r_circ_opt);
fprintf('Delta V2       : %.6f km/s (%.2f m/s)\n', dv_circ_min, dv_circ_min * 1000);

fprintf('\n>>> 2. MIGLIOR ORBITA ELLITTICA\n');
fprintf('Raggio Pericentro : %.4f km\n', r_p_ell_opt);
fprintf('Raggio Apocentro  : %.4f km\n', r_a_ell_opt);
fprintf('Semiasse Maggiore : %.4f km\n', a_ell_opt);
fprintf('Eccentricità      : %.4f\n', e_ell_opt);
fprintf('Delta V2          : %.6f km/s (%.2f m/s)\n', dv_ell_min, dv_ell_min * 1000);

fprintf('\n>>> 3. VERDETTO GLOBALE DELLA MISSIONE\n');
if dv_ell_min < dv_circ_min
    fprintf('La strategia migliore in assoluto è: ORBITA ELLITTICA\n');
    DV2_best = dv_ell_min;
    r_p_best = r_p_ell_opt;
    e_h_best = e_h_ell_opt;
else
    fprintf('La strategia migliore in assoluto è: ORBITA CIRCOLARE\n');
    DV2_best = dv_circ_min;
    r_p_best = r_circ_opt;
    e_h_best = e_h_circ_opt;
end

% Calcolo Delta V Totale (presumendo che DeltaV_1 sia nel workspace)
if exist('DeltaV_1', 'var')
    DV_tot = DeltaV_1 + DV2_best;
    fprintf('\nDELTA V1 (Partenza Terra) : %.4f km/s\n', DeltaV_1);
    fprintf('DELTA V2 (Arrivo Ast.)    : %.4f km/s\n', DV2_best);
    fprintf('------------------------------------------------------\n');
    fprintf('DELTA V TOTALE            : %.4f km/s\n', DV_tot);
else
    fprintf('\nNota: DeltaV_1 non trovato nel workspace per calcolare il Totale.\n');
end
fprintf('======================================================\n');


% =========================================================================
% --- PLOT GRAFICI DEL DELTA V ---
% =========================================================================

% 1. Plot DV vs Raggio per orbite circolari
figure('Name', 'Analisi Costo Orbita Circolare');
plot(orb_circ.r_ph, orb_circ.dv * 1000, 'b', 'LineWidth', 2);
grid on;
title('Costo Cattura in Orbita Circolare');
xlabel('Raggio dell''orbita [km]');
ylabel('\DeltaV [m/s]');

% 2. Plot 3D (Surface) DV vs Eccentricità e Semiasse per orbite ellittiche
figure('Name', 'Surface Costo Orbita Ellittica');
surf(orb_ell.e, orb_ell.a, orb_ell.dv * 1000, 'EdgeColor', 'none');
colorbar;
grid on; view(3);
title('Costo Cattura in Orbita Ellittica (Surface)');
xlabel('Eccentricità [-]');
ylabel('Semiasse Maggiore a [km]');
zlabel('\DeltaV [m/s]');

% 3. Plot 2D DV vs Semiasse a parità di Eccentricità Ottima
figure('Name', 'Analisi Costo a Eccentricità Fissa');
% Estraggo tutta la colonna corrispondente all'eccentricità ottima
dv_slice = orb_ell.dv(:, col_ell);
a_slice  = orb_ell.a(:, col_ell);
plot(a_slice, dv_slice * 1000, 'r', 'LineWidth', 2);
grid on;
title(sprintf('Costo Cattura Ellittica a Eccentricità fissa (e = %.4f)', e_ell_opt));
xlabel('Semiasse Maggiore a [km]');
ylabel('\DeltaV [m/s]');


% =========================================================================
% --- PLOT ORBITE NEL SISTEMA CENTRATO SULL'ASTEROIDE (3D e 2D) ---
% =========================================================================

% Preparazione vettori angolo
th_ell = linspace(0, 2*pi, 200);
th_inf = acos(-1/e_h_best);
th_hyp = linspace(-th_inf + 0.05, th_inf - 0.05, 200);

% Parametri Iperbole d'arrivo (basata sulla manovra migliore)
r_hyp_mag = a_H2 * (1 - e_h_best^2) ./ (1 + e_h_best * cos(th_hyp));

% Parametri Orbita Circolare Ottima
r_circ_mag = ones(size(th_ell)) * r_circ_opt;

% Parametri Orbita Ellittica Ottima
p_ell = a_ell_opt * (1 - e_ell_opt^2);
r_ell_mag = p_ell ./ (1 + e_ell_opt * cos(th_ell));

% --- PLOT 3D ---
figure('Name', 'Scenario 4: Arrivo sull''Asteroide (3D)');
hold on; grid on; axis equal; view(3);

% Disegno Asteroide
[X_A, Y_A, Z_A] = sphere(40);
surf(X_A*(diam_ast/2), Y_A*(diam_ast/2), Z_A*(diam_ast/2), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'DisplayName', 'Asteroide 2003 UC20');

% Iperbole Arrivo
plot3(r_hyp_mag .* cos(th_hyp), r_hyp_mag .* sin(th_hyp), zeros(size(th_hyp)), 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Arrivo');
% Circolare
plot3(r_circ_mag .* cos(th_ell), r_circ_mag .* sin(th_ell), zeros(size(th_ell)), '--b', 'LineWidth', 1.5, 'DisplayName', 'Orbita Circ. Ottimale');
% Ellittica
plot3(r_ell_mag .* cos(th_ell), r_ell_mag .* sin(th_ell), zeros(size(th_ell)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Orbita Ellittica Ottimale');

% Marker Punto di Manovra (Pericentro)
plot3(r_p_best, 0, 0, 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV: Punto di Cattura');

lim_zoom_ast = max(r_a_ell_opt, r_circ_opt) * 1.5;
xlim([-lim_zoom_ast lim_zoom_ast]); ylim([-lim_zoom_ast lim_zoom_ast]); zlim([-lim_zoom_ast lim_zoom_ast]);
xlabel('Asse Pericentrale P [km]'); ylabel('Asse Trasverso Q [km]'); zlabel('Z [km]');
title('Sistema Centrato sull''Asteroide (3D)');
legend('Location', 'best');
hold off;

% --- PLOT 2D ---
figure('Name', 'Scenario 4: Arrivo sull''Asteroide (2D)');
hold on; grid on; axis equal; view(2);

% Disegno Asteroide in 2D (come cerchio)
fill((diam_ast/2)*cos(th_ell), (diam_ast/2)*sin(th_ell), [0.5 0.5 0.5], 'DisplayName', 'Asteroide');

plot(r_hyp_mag .* cos(th_hyp), r_hyp_mag .* sin(th_hyp), 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Arrivo');
plot(r_circ_mag .* cos(th_ell), r_circ_mag .* sin(th_ell), '--b', 'LineWidth', 1.5, 'DisplayName', 'Orbita Circ. Ottimale');
plot(r_ell_mag .* cos(th_ell), r_ell_mag .* sin(th_ell), 'g', 'LineWidth', 1.5, 'DisplayName', 'Orbita Ellittica Ottimale');
plot(r_p_best, 0, 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV: Punto di Cattura');

xlim([-lim_zoom_ast lim_zoom_ast]); ylim([-lim_zoom_ast lim_zoom_ast]);
xlabel('Asse Pericentrale P [km]'); ylabel('Asse Trasverso Q [km]');
title('Vista 2D Piano Orbitale Asteroide');
legend('Location', 'best');
hold off;