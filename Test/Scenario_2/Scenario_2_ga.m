clear
close all
clc

% --- INIZIALIZZAZIONE STORICO ---
global history_ga history_fmincon;
history_ga = [];
history_fmincon = [];

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
% --- FASE 1: RICERCA GLOBALE CON GA --- 
fprintf('Fase 1: Avvio Algoritmo Genetico...\n');

ga_opts = optimoptions('ga', 'PopulationSize', 200, 'MaxGenerations', 100, 'Display', 'iter', 'OutputFcn', @outfun_ga);

% Definiamo la funzione anonima per passare i dati dell'asteroide
obj_fun = @(x) objective_function(x, ast);
% Non passo @(x)constraints per alleggerire il codice, altrimenti si
% dovrebbe calcolare una funzione di vincolo per ogni individuo
[x_ga, fval_ga] = ga(obj_fun, 3, [], [], [], [], lb, ub, [], ga_opts);
% --- FASE 2: RIFINITURA LOCALE CON FMINCON --- 
fprintf('\nFase 2: Rifinitura con fmincon (Partendo dal risultato GA)...\n');

fm_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'TolFun', 1e-8, 'OutputFcn', @outfun_fmincon);

% Usiamo x_ga come punto di partenza (x0)
[x_opt, dv_opt] = fmincon(obj_fun, x_ga, [], [], [], [], lb, ub, @(x) constraints(x, ast), fm_opts);
% 3. RISULTATI FINALI 
% Stampa dei 3 parametri ottimi trovati
fprintf('\n--- RISULTATI OTTIMI ---');
fprintf('\nDelta V Totale: %.4f km/s', dv_opt);
fprintf('\nRiga 1: Theta_1 Partenza (Terra)        = %.2f deg', rad2deg(x_opt(1)));
fprintf('\nRiga 2: Theta_2 Arrivo (Asteroide)      = %.2f deg', rad2deg(x_opt(2)));
fprintf('\nRiga 3: Omega_T Argomento Pericentro    = %.2f deg\n', rad2deg(x_opt(3)));

% =========================================================================
% 4. PLOTTING DELLE CURVE DI CONVERGENZA
% =========================================================================
fprintf('\nGenerazione dei grafici di convergenza...\n');

% --- FIGURE 1: Convergenza Algoritmo Genetico (GA) ---
if ~isempty(history_ga)
    figure('Name', 'Convergenza GA', 'Color', 'w', 'Position', [100, 100, 700, 500]);
    
    pop_size = ga_opts.PopulationSize; 
    num_righe = size(history_ga, 1);
    num_gen = floor(num_righe / pop_size);
    
    best_dv_ga = zeros(num_gen, 1);
    mean_dv_ga = zeros(num_gen, 1);
    all_dv_ga = zeros(num_gen, pop_size); % Matrice per tracciare tutte le curve
    
    for g = 1:num_gen
        idx_start = (g-1)*pop_size + 1;
        idx_end = g*pop_size;
        
        scores_gen = history_ga(idx_start:idx_end, 4);
        
        % Ordiniamo i punteggi della generazione per formare curve di esplorazione
        scores_sorted = sort(scores_gen);
        all_dv_ga(g, :) = scores_sorted;
        
        best_dv_ga(g) = scores_sorted(1); 
        
        valid_scores = scores_gen(scores_gen < 1e5);
        if isempty(valid_scores)
            mean_dv_ga(g) = NaN;
        else
            mean_dv_ga(g) = mean(valid_scores);
        end
    end
    
    hold on;
    
    % Sostituiamo i valori scartati per penalità (es. 1e6) con NaN in modo che 
    % le linee grigie non distruggano la scala dell'asse Y del grafico
    all_dv_ga(all_dv_ga >= 1e5) = NaN;
    
    % Plot di tutti gli individui in grigio chiaro (sottofondo)
    h_altri = plot(1:num_gen, all_dv_ga, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    
    % Plot del migliore e della media sopra le linee grigie (Media ora continua '-c')
    h_best = plot(1:num_gen, best_dv_ga, '-b', 'LineWidth', 2.5);
    h_mean = plot(1:num_gen, mean_dv_ga, '-c', 'LineWidth', 2);
    
    [min_val, min_idx] = min(best_dv_ga);
    % h_opt = plot(min_idx, min_val, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
    
    grid on;
    title('Fase 1: Convergenza Algoritmo Genetico');
    xlabel('Generazioni');
    ylabel('\DeltaV Totale [km/s]');
    
    % Legenda personalizzata per raggruppare tutte le linee grigie sotto una sola voce
    legend([h_best, h_mean, h_altri(1)], ...
        {'Miglior Individuo', 'Media Popolazione', 'Altri Individui'}, 'Location', 'northeast');
end

% --- FIGURE 2: Convergenza fmincon ---
if ~isempty(history_fmincon)
    figure('Name', 'Convergenza fmincon', 'Color', 'w', 'Position', [150, 150, 700, 500]);
    
    iters_fmincon = 1:size(history_fmincon, 1);
    dv_fmincon = history_fmincon(:, 4);
    
    plot(iters_fmincon, dv_fmincon, '-o', 'Color', '#D95319', 'LineWidth', 1.5, ...
        'MarkerFaceColor', '#EDB120', 'DisplayName', 'Evoluzione \DeltaV');
    hold on;
    
    plot(iters_fmincon(1), dv_fmincon(1), 'sg', 'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
        'DisplayName', 'Partenza (Risultato GA)');
    plot(iters_fmincon(end), dv_fmincon(end), 'pr', 'MarkerSize', 14, 'MarkerFaceColor', 'y', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Ottimo Locale Finale');
    
    grid on;
    title('Fase 2: Rifinitura fmincon');
    xlabel('Iterazioni');
    ylabel('\DeltaV Totale [km/s]');
    legend('Location', 'northeast');
    
    if max(dv_fmincon) - min(dv_fmincon) < 1e-4
        margine = 1e-4;
        ylim([min(dv_fmincon)-margine, max(dv_fmincon)+margine]);
    end
end