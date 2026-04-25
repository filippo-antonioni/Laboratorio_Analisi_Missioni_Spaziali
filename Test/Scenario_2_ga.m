clear
close all
clc

mu = 398600;

% Dati dell'asteroide 363505 (2003 UC20) 

% Costante di conversione
AU_to_km = 149597870.7; % Unità Astronomica in km

% --- Parametri Fisici ---
M = 6.9140e12;            % Massa in [kg]
D = 1.88;                 % Diametro in [km]

% --- Parametri Orbitali ---
AU = 149597870.7;
ast.a = 0.781241 * AU; 
ast.e = 0.336932;
ast.i = deg2rad(3.78);
ast.OM = deg2rad(187.92);
ast.om = deg2rad(60.16);

% 2. CONFIGURAZIONE OTTIMIZZAZIONE
% Variabili x = [theta1_partenza, theta2_arrivo, omega_trasferimento]
lb = [0, 0, 0];
ub = [2*pi, 2*pi, 2*pi];

% --- FASE 1: RICERCA GLOBALE CON GA --- [cite: 120]
fprintf('Fase 1: Avvio Algoritmo Genetico...\n');
ga_opts = optimoptions('ga', 'PopulationSize', 200, 'MaxGenerations', 100, 'Display', 'iter');

% Definiamo la funzione anonima per passare i dati dell'asteroide
obj_fun = @(x) objective_function(x, ast);

% Non passo @(x)constraints per alleggerire il codice, altrimenti si
% dovrebbe calcolare una funzione di vincolo per ogni individuo
[x_ga, fval_ga] = ga(obj_fun, 3, [], [], [], [], lb, ub, [], ga_opts);

% --- FASE 2: RIFINITURA LOCALE CON FMINCON --- 
fprintf('\nFase 2: Rifinitura con fmincon (Partendo dal risultato GA)...\n');
fm_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'TolFun', 1e-8);

% Usiamo x_ga come punto di partenza (x0)
[x_opt, dv_opt] = fmincon(obj_fun, x_ga, [], [], [], [], lb, ub, @(x) constraints(x, ast), fm_opts);

% 3. RISULTATI FINALI 
% Stampa dei 3 parametri ottimi trovati
fprintf('\n--- RISULTATI OTTIMI ---');
fprintf('\nDelta V Totale: %.4f km/s', dv_opt);
fprintf('\nRiga 1: Theta_1 Partenza (Terra)        = %.2f deg', rad2deg(x_opt(1)));
fprintf('\nRiga 2: Theta_2 Arrivo (Asteroide)      = %.2f deg', rad2deg(x_opt(2)));
fprintf('\nRiga 3: Omega_T Argomento Pericentro    = %.2f deg\n', rad2deg(x_opt(3)));









