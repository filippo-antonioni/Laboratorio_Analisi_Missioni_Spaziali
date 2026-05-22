clear
clc
close all

%% SECTION 1 - Geocentrica -> Eliocentrica (Iperbole di Fuga)
%%% Load dati orbita finale geocentrica (scenario 1) - ORBITA DI PARCHEGGIO %%%
T = load("DatiSC1-2026.txt");
mu_earth = 398600;
m_earth = 5.974e24; % [kg]
gruppo = 23;
R_earth = 6371;

% Orbita target/finale (vettori di stato)
rr_f = T(gruppo,8:10)';
vv_f = T(gruppo,11:13)';
state_f = [rr_f; vv_f];

% Trasformazione per avere i parametri kepleriani bersaglio
[a_op, e_op, i_op, OM_op, om_op, th_op] = car2par(rr_f, vv_f, mu_earth);

% Moduli del raggio e velocità di pericentro dell'orbita di parcheggio
r_op_p = a_op * (1 - e_op);
v_op_p = sqrt(mu_earth * (2/r_op_p - 1/a_op));
r_op_a = a_op * (1 + e_op);

%%% Load dati orbita eliocentrica (scenario 2) %%%
load('workspace_sc2_gridsearch.mat'); 
mu_sun  = 1.32712440018e11; % Parametro gravitazionale del Sole [km^3/s^2]
m_sun   = 1.989e30;         % [kg]
R_sun   = 696340;           % Raggio del Sole [km]
AU      = 149597870.7;      % 1 Unità Astronomica [km]

%%% Calcolo SOI terra (usando la distanza reale alla partenza) %%%
R_terra_partenza = norm(r1_opt); 
r_SOI_terra = R_terra_partenza * (m_earth/m_sun)^(2/5);

% Velocità eliocentrica della Terra alla partenza (V1i)
[r1i, V1i] = par2car(a_T, e_T, i_T, OM_T, om_T, ottimo_th1, mu_sun);

% Velocità eliocentrica del satellite sull'orbita di trasferimento alla partenza (V1T)
[R1T, V1T] = par2car(ottimo_aT, ottimo_eT, i_trasf_opt, OM_transf_opt, ottimo_omT, th1_T_opt, mu_sun);

% Eccesso iperbolico vettoriale e modulo
v_inf_1_vec = V1T - V1i;

% Matrice di rotazione da SDR eclittico a SDR ECI
epsilon = deg2rad(23.45);
T_rot = [1 0 0;
         0 cos(epsilon) sin(epsilon);
         0 -sin(epsilon) cos(epsilon)];
         
% Converto la v_inf da eclittico a ECI 
v_inf_1_vec = T_rot' * v_inf_1_vec;
v_inf_1 = norm(v_inf_1_vec);

% Trovo versore che definisce l'asintoto uscente dell'iperbole 
r_inf_vers = v_inf_1_vec / v_inf_1;

% Procedo a calcolare con un for il raggio r_h
th_h_op_vec = linspace(0, 2*pi, 1000);
N_r_1 = length(th_h_op_vec); 

% Calcolo semiasse maggiore iperbole
a_h = -mu_earth / v_inf_1^2;

% Vettori che uso nel ciclo for
dv_1_vec = zeros(N_r_1, 1);
r_h_vec = zeros(N_r_1, 3);
e_h_vec = zeros(N_r_1, 1);
th_h_vec = zeros(N_r_1, 1);
th_inf_vec = zeros(N_r_1, 1);

kk = [0 0 1]';
margine_terra = 200;

% Opzioni per fsolve 
options = optimoptions('fsolve', 'Display', 'none');

for i = 1:N_r_1 
    % Calcolo r_h_i associato al th_h_i della parking orbit
    [r_h_i, v_h_i] = par2car(a_op, e_op, i_op, OM_op, om_op, th_h_op_vec(i), mu_earth);
    
    r_h_i_norm = norm(r_h_i);
    r_h_i_vers = r_h_i / r_h_i_norm;
    alfa = acos(dot(r_h_i_vers, r_inf_vers));
    
    % Calcolo e_h guess iniziale
    e_h_0 = 1 - r_h_i_norm/a_h;
    th_inf_0 = acos(-1 / e_h_0);
    th_h_0 = th_inf_0 - alfa;
    
    x0 = [e_h_0; th_h_0; th_inf_0];
    
    fun = @(x) [r_h_i_norm - a_h*(1-x(1)^2)/(1+x(1)*cos(x(2)));
                cos(x(3)) + 1/x(1);
                x(3) - x(2) - alfa];   
    
    [x_i, ~, exitflag] = fsolve(fun, x0, options); 
    
    % Calcolo raggio di pericentro iperbole i-esima
    r_p_h_i = a_h*(1 - x_i(1));
    
    % Controllo vincoli
    if exitflag > 0 && (x_i(3) > pi/2 && x_i(3) < pi) && r_p_h_i > (R_earth + margine_terra)
        
        r_h_vec(i,:) = r_h_i;
        e_h_vec(i) = x_i(1);
        th_h_vec(i) = x_i(2);
        th_inf_vec(i) = x_i(3); 
        
        % Direzione momento angolare specifico
        h_i_vec = cross(r_h_i, r_inf_vers);
        h_i_vers = h_i_vec / norm(h_i_vec);
        
        % Calcolo inclinazione e RAAN
        i_h_i = acos(h_i_vers(3));
        NN_h_i = cross(kk, h_i_vers) / norm(cross(kk, h_i_vers));
        if NN_h_i(2) >= 0
            OM_h_i = acos(NN_h_i(1)); 
        else 
            OM_h_i = 2*pi - acos(NN_h_i(1)); 
        end
        
        cos_beta_i = dot(NN_h_i, r_h_i_vers);
        sin_beta_i = dot(cross(NN_h_i, r_h_i_vers), h_i_vers);
        beta_i = atan2(sin_beta_i, cos_beta_i);
        
        om_h_i = beta_i - x_i(2); 
        
        % calcolo parametri (vettori di stato sull'iperbole in quel punto)
        [~, vv_h_i] = par2car(a_h, x_i(1), i_h_i, OM_h_i, om_h_i, x_i(2), mu_earth);
        
        dv_1_vec(i) = norm(vv_h_i - v_h_i);
    else 
        % Caso in cui non rispetto i vincoli
        r_h_vec(i,:) = r_h_i;
        e_h_vec(i) = inf;
        th_h_vec(i) = inf;
        th_inf_vec(i) = inf; 
        dv_1_vec(i) = inf;
    end 
end

% Trovo il minimo 
[dv_opt_1, idx_opt] = min(dv_1_vec);

%% --- RECUPERO PARAMETRI OTTIMI PER I PLOT ---
e_H_opt = e_h_vec(idx_opt);
th_H_opt = th_h_vec(idx_opt);
th_inf_opt = th_inf_vec(idx_opt);
r_manovra_opt = r_h_vec(idx_opt, :)';

% Ricalcolo gli angoli orbitali per l'iperbole ottima (necessari per il plot 3D)
h_opt_vec = cross(r_manovra_opt, r_inf_vers);
h_opt_vers = h_opt_vec / norm(h_opt_vec);
i_H_opt = acos(h_opt_vers(3));

NN_H_opt = cross(kk, h_opt_vers) / norm(cross(kk, h_opt_vers));
if NN_H_opt(2) >= 0
    OM_H_opt = acos(NN_H_opt(1)); 
else 
    OM_H_opt = 2*pi - acos(NN_H_opt(1)); 
end

cos_beta_opt = dot(NN_H_opt, r_manovra_opt / norm(r_manovra_opt));
sin_beta_opt = dot(cross(NN_H_opt, r_manovra_opt / norm(r_manovra_opt)), h_opt_vers);
beta_opt = atan2(sin_beta_opt, cos_beta_opt);
om_H_opt = beta_opt - th_H_opt;

%%% 4. CONTROLLO DI COERENZA DELL'IPERBOLE %%%
if e_H_opt > 1 && a_h < 0
    fprintf('\n[OK] Controllo superato: L''orbita di fuga è un''iperbole (e = %.4f, a = %.2f km).\n', e_H_opt, a_h);
    fprintf('Costo manovra ottimale (Delta V1): %.4f km/s\n', dv_opt_1);
else
    warning('ATTENZIONE: I parametri calcolati non corrispondono a un''iperbole!');
end

%%% 5. PLOT DEL SISTEMA GEOCENTRICO 3D %%%
figure('Name', 'Scenario 3: Iperbole di Fuga dalla Terra (3D)');
hold on; grid on; axis equal; view(3);

[X_E, Y_E, Z_E] = sphere(50);
load topo; 
surf(X_E*R_earth, Y_E*R_earth, Z_E*R_earth, 'CData', topo, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', 'HandleVisibility', 'off');
colormap(topomap1); 
plot3(NaN, NaN, NaN, 'o', 'Color', 'b', 'MarkerFaceColor', [0.2 0.5 0.8], ...
    'MarkerSize', 8, 'DisplayName', 'Terra');

th_ell = linspace(0, 2*pi, 200);
th_hyp = linspace(-th_inf_opt + 0.05, th_inf_opt - 0.05, 200);

r_op_plot = zeros(3, length(th_ell));
for k = 1:length(th_ell)
    [r_op_plot(:,k), ~] = par2car(a_op, e_op, i_op, OM_op, om_op, th_ell(k), mu_earth);
end
plot3(r_op_plot(1,:), r_op_plot(2,:), r_op_plot(3,:), 'b', 'LineWidth', 1.5, 'DisplayName', 'Orbita di Parcheggio');

r_hyp_plot = zeros(3, length(th_hyp));
for k = 1:length(th_hyp)
    [r_hyp_plot(:,k), ~] = par2car(a_h, e_H_opt, i_H_opt, OM_H_opt, om_H_opt, th_hyp(k), mu_earth);
end
plot3(r_hyp_plot(1,:), r_hyp_plot(2,:), r_hyp_plot(3,:), 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Fuga');

plot3(r_manovra_opt(1), r_manovra_opt(2), r_manovra_opt(3), 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV1: Punto di Manovra');

lim_zoom = r_op_a * 2; 
xlim([-lim_zoom lim_zoom]); ylim([-lim_zoom lim_zoom]); zlim([-lim_zoom lim_zoom]);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Iniezione su Iperbole di Fuga Geocentrica (Vista 3D)');
legend('Location', 'best');
hold off;

%%% 6. PLOT DEL SISTEMA GEOCENTRICO 2D %%%
figure('Name', 'Scenario 3: Iperbole di Fuga (2D Dimensione e Forma)');
hold on; grid on; axis equal;
[X_E2, Y_E2, Z_E2] = sphere(50);
surf(X_E2*R_earth, Y_E2*R_earth, Z_E2*R_earth, 'CData', topo, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', 'HandleVisibility', 'off');
colormap(topomap1);

plot3(NaN, NaN, NaN, 'o', 'Color', 'b', 'MarkerFaceColor', [0.2 0.5 0.8], ...
    'MarkerSize', 8, 'DisplayName', 'Terra');

r_op_mag = a_op * (1 - e_op^2) ./ (1 + e_op * cos(th_ell));
z_layer_ell = ones(size(th_ell)) * R_earth; 
plot3(r_op_mag .* cos(th_ell), r_op_mag .* sin(th_ell), z_layer_ell, 'b', 'LineWidth', 1.5, 'DisplayName', 'Orbita di Parcheggio');

r_hyp_mag = a_h * (1 - e_H_opt^2) ./ (1 + e_H_opt * cos(th_hyp));
z_layer_hyp = ones(size(th_hyp)) * R_earth;
plot3(r_hyp_mag .* cos(th_hyp), r_hyp_mag .* sin(th_hyp), z_layer_hyp, 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Fuga');

r_manovra_mag = norm(r_manovra_opt);
plot3(r_manovra_mag*cos(th_H_opt), r_manovra_mag*sin(th_H_opt), R_earth, 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV1: Punto di Manovra');

view(2);
xlim([-lim_zoom lim_zoom]); ylim([-lim_zoom lim_zoom]);
xlabel('Asse Pericentrale P [km]'); ylabel('Asse Trasverso Q [km]');
title('Vista 2D Sagome Geometriche (Terra)');
legend('Location', 'best');
hold off;

%%% 7. CALCOLO TOF %%%
cos_th_SOI_1 = (a_h * (1 - e_H_opt^2) / r_SOI_terra - 1) / e_H_opt;
cos_th_SOI_1 = max(-1, min(1, cos_th_SOI_1)); 

H_1 = acosh((e_H_opt + cos_th_SOI_1) / (1 + e_H_opt * cos_th_SOI_1));
TOF_fuga_sec = sqrt((-a_h)^3 / mu_earth) * (e_H_opt * sinh(H_1) - H_1);
TOF_fuga_giorni = TOF_fuga_sec / (24 * 3600);
fprintf('\nTOF Fuga dalla Terra (da Manovra a SOI): %.2f giorni (%.0f sec)\n', TOF_fuga_giorni, TOF_fuga_sec);


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
mu_ast = G * m_ast / (1e9); % diviso per 1e9 per portare in km^3
a_H2 = - mu_ast / v_inf_2^2;

%%% CONTROLLO SEMIASSE MAGGIORE DELL'IPERBOLE %%%
if a_H2 < 0
    fprintf('\n[OK] Controllo superato: L''orbita di arrivo è un''iperbole (a = %.2e km).\n', a_H2);
else
    error('ATTENZIONE: I parametri calcolati non corrispondono a un''iperbole!');
end

% Creazione vettori di scansione
margine = 0.4; % Aggiunto punto e virgola mancante
rvec = linspace(diam_ast/2 + margine, r_SOI_ast, 100);
e_vec = linspace(0, 1, 1000);
e_vec = e_vec(2:end-1); % Rimuovo 0 (è il caso circolare) e 1 (parabola)

N_r = length(rvec);
N_e = length(e_vec);

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
    r_a_k = r_k;
    
    for j = 1:N_e
        e_k = e_vec(j);
        
        % Calcolo subito il pericentro
        r_p_k = r_a_k * (1 - e_k) / (1 + e_k);
        
        % check per evitare che orbita mi attraversi l'asteroide
        [orbita_ok, ~] = check_feasibility(r_p_k, r_a_k, diam_ast/2, margine);
        
        if orbita_ok
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
                orb_ell.a(i,j)      = a_k;
                orb_ell.e_h(i,j)    = e_h_ell;
                orb_ell.dv(i, j)    = abs(v_p_k - v_ph_ell);
                orb_ell.delta(i, j) = -a_H2 * sqrt(e_h_ell^2 - 1);
            else
                orb_ell.dv(i, j) = NaN;
            end
        else
            % Scarto: L'orbita NON ha passato il controllo geometrico (Impatto, Fuga, ecc.)
            orb_ell.dv(i, j) = NaN;
        end
    end
end

% =========================================================================
% --- ESTRAZIONE DEI RISULTATI OTTIMALI ---
% =========================================================================
% 1. Ottimo per Orbita Circolare
[dv_circ_min, idx_circ] = min(orb_circ.dv);
r_circ_opt = orb_circ.r_ph(idx_circ);
e_h_circ_opt = orb_circ.e_h(idx_circ);

% 2. Ottimo per Orbita Ellittica
[dv_ell_min, idx_ell_linear] = min(orb_ell.dv(:));
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
    fprintf('La strategia migliore in assoluto all''arrivo è: ORBITA ELLITTICA\n');
    DV2_best = dv_ell_min;
    r_p_best = r_p_ell_opt;
    e_h_best = e_h_ell_opt;
else
    fprintf('La strategia migliore in assoluto all''arrivo è: ORBITA CIRCOLARE\n');
    DV2_best = dv_circ_min;
    r_p_best = r_circ_opt;
    e_h_best = e_h_circ_opt;
end

% Calcolo Delta V Totale (Modificato con dv_opt_1 calcolato nella Sezione 1)
if exist('dv_opt_1', 'var')
    DV_tot = dv_opt_1 + DV2_best;
    fprintf('\nDELTA V1 (Partenza Terra) : %.4f km/s\n', dv_opt_1);
    fprintf('DELTA V2 (Arrivo Ast.)    : %.4f km/s\n', DV2_best);
    fprintf('------------------------------------------------------\n');
    fprintf('DELTA V TOTALE DELLA MISS.: %.4f km/s\n', DV_tot);
else
    fprintf('\nNota: dv_opt_1 non trovato nel workspace per calcolare il Totale.\n');
end
fprintf('======================================================\n');

% =========================================================================
% --- PLOT GRAFICI DEL DELTA V ---
% =========================================================================
% 1. Plot DV vs Raggio per orbite circolari
figure('Name', 'Analisi Costo Orbita Circolare');
plot(orb_circ.r_ph, orb_circ.dv * 1000, '-bo', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
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
dv_slice = orb_ell.dv(:, col_ell);
a_slice  = orb_ell.a(:, col_ell);
valid_idx = ~isnan(dv_slice);
plot(a_slice(valid_idx), dv_slice(valid_idx) * 1000, '-rs', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
grid on;
title(sprintf('Costo Cattura Ellittica a Eccentricità fissa (e = %.4f)', e_ell_opt));
xlabel('Semiasse Maggiore a [km]');
ylabel('\DeltaV [m/s]');

% =========================================================================
% --- PLOT ORBITE NEL SISTEMA CENTRATO SULL'ASTEROIDE (3D e 2D) ---
% =========================================================================
th_ell = linspace(0, 2*pi, 200);
r_limite_plot = r_SOI_ast * 2; 
arg_cos = (a_H2 * (1 - e_h_best^2) / r_limite_plot - 1) / e_h_best;
arg_cos = max(-1, min(1, arg_cos)); 
th_limite = acos(arg_cos);

th_hyp_ast = linspace(-th_limite, th_limite, 300);

r_hyp_mag_ast = a_H2 * (1 - e_h_best^2) ./ (1 + e_h_best * cos(th_hyp_ast));
r_circ_mag = ones(size(th_ell)) * r_circ_opt;
p_ell = a_ell_opt * (1 - e_ell_opt^2);
r_ell_mag = p_ell ./ (1 + e_ell_opt * cos(th_ell));

% --- PLOT 3D ---
figure('Name', 'Scenario 4: Arrivo sull''Asteroide (3D)');
hold on; grid on; axis equal; view(3);
[X_A, Y_A, Z_A] = sphere(40);
surf(X_A*(diam_ast/2), Y_A*(diam_ast/2), Z_A*(diam_ast/2), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'DisplayName', 'Asteroide 2003 UC20');

plot3(r_hyp_mag_ast .* cos(th_hyp_ast), r_hyp_mag_ast .* sin(th_hyp_ast), zeros(size(th_hyp_ast)), 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Arrivo');
plot3(r_circ_mag .* cos(th_ell), r_circ_mag .* sin(th_ell), zeros(size(th_ell)), '--b', 'LineWidth', 1.5, 'DisplayName', 'Orbita Circ. Ottimale');
plot3(r_ell_mag .* cos(th_ell), r_ell_mag .* sin(th_ell), zeros(size(th_ell)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Orbita Ellittica Ottimale');
plot3(r_p_best, 0, 0, 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV2: Punto di Cattura');

lim_zoom_ast = max(r_a_ell_opt, r_circ_opt) * 1.5;
xlim([-lim_zoom_ast lim_zoom_ast]); ylim([-lim_zoom_ast lim_zoom_ast]); zlim([-lim_zoom_ast lim_zoom_ast]);
xlabel('Asse Pericentrale P [km]'); ylabel('Asse Trasverso Q [km]'); zlabel('Z [km]');
title('Sistema Centrato sull''Asteroide (3D)');
legend('Location', 'best');
hold off;

% --- PLOT 2D ---
figure('Name', 'Scenario 4: Arrivo sull''Asteroide (2D)');
hold on; grid on; axis equal; view(2);
fill((diam_ast/2)*cos(th_ell), (diam_ast/2)*sin(th_ell), [0.5 0.5 0.5], 'DisplayName', 'Asteroide');

plot(r_hyp_mag_ast .* cos(th_hyp_ast), r_hyp_mag_ast .* sin(th_hyp_ast), 'r', 'LineWidth', 2, 'DisplayName', 'Iperbole di Arrivo');
plot(r_circ_mag .* cos(th_ell), r_circ_mag .* sin(th_ell), '--b', 'LineWidth', 1.5, 'DisplayName', 'Orbita Circ. Ottimale');
plot(r_ell_mag .* cos(th_ell), r_ell_mag .* sin(th_ell), 'g', 'LineWidth', 1.5, 'DisplayName', 'Orbita Ellittica Ottimale');
plot(r_p_best, 0, 'pk', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', '\DeltaV2: Punto di Cattura');

xlim([-lim_zoom_ast lim_zoom_ast]); ylim([-lim_zoom_ast lim_zoom_ast]);
xlabel('Asse Pericentrale P [km]'); ylabel('Asse Trasverso Q [km]');
title('Vista 2D Piano Orbitale Asteroide');
legend('Location', 'best');
hold off;

%% TOF PER IPERBOLE DI ARRIVO
cos_th_SOI_2 = (a_H2 * (1 - e_h_best^2) / r_SOI_ast - 1) / e_h_best;
cos_th_SOI_2 = max(-1, min(1, cos_th_SOI_2));

H_2 = acosh((e_h_best + cos_th_SOI_2) / (1 + e_h_best * cos_th_SOI_2));
TOF_arrivo_sec = sqrt((-a_H2)^3 / mu_ast) * (e_h_best * sinh(H_2) - H_2);
TOF_arrivo_giorni = TOF_arrivo_sec / (24 * 3600);
fprintf('\nTOF Arrivo su Asteroide (da SOI a Pericentro): %.2f giorni (%.0f sec)\n', TOF_arrivo_giorni, TOF_arrivo_sec);

%% TEMPO DELLA MANOVRA CIRCOLARE/ELLITTICA (NON OTTIMALE)
fprintf('\n--- CONFRONTO TEMPI DI VOLO (STRATEGIA SCARTATA) ---\n');
if dv_ell_min < dv_circ_min
    e_h_loser = e_h_circ_opt;
    nome_loser = 'CIRCOLARE';
else
    e_h_loser = e_h_ell_opt;
    nome_loser = 'ELLITTICA';
end

cos_th_SOI_loser = (a_H2 * (1 - e_h_loser^2) / r_SOI_ast - 1) / e_h_loser;
cos_th_SOI_loser = max(-1, min(1, cos_th_SOI_loser)); 

H_loser = acosh((e_h_loser + cos_th_SOI_loser) / (1 + e_h_loser * cos_th_SOI_loser));
TOF_loser_sec = sqrt((-a_H2)^3 / mu_ast) * (e_h_loser * sinh(H_loser) - H_loser);
TOF_loser_giorni = TOF_loser_sec / (24 * 3600);
fprintf('TOF Arrivo su Asteroide per orbita %s (Scartata): %.2f giorni (%.0f sec)\n', nome_loser, TOF_loser_giorni, TOF_loser_sec);
fprintf('======================================================\n');