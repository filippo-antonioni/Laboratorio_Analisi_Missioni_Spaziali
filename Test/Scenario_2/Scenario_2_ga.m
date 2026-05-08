clear
close all
clc

% --- INIZIALIZZAZIONE STORICO ---
global history_ga history_fmincon;
history_ga = [];
history_fmincon = [];

best_dv_global = inf;
best_x_global = [];
best_history_ga = [];
best_history_fmincon = [];

N_runs = 200;
results_x = zeros(N_runs, 3);
results_dv = zeros(N_runs, 1);
results_time = zeros(N_runs, 1);

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

%Inizio ciclo for della run
for k=1:N_runs
   
    %azzero variabili per algoritmi genetici
    history_ga=[];
    history_fmincon=[];
    tic %inizio a contare il tempo
    
    % --- FASE 1: RICERCA GLOBALE CON GA --- 
    ga_opts = optimoptions('ga', 'PopulationSize', 200, 'MaxGenerations', 100, 'Display', 'off', 'OutputFcn', @outfun_ga);
    obj_fun = @(x) objective_function(x, ast);
    [x_ga, fval_ga] = ga(obj_fun, 3, [], [], [], [], lb, ub, [], ga_opts);
    
    % --- FASE 2: RIFINITURA LOCALE CON FMINCON --- 
    fm_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'TolFun', 1e-8, 'OutputFcn', @outfun_fmincon);
    [x_opt, dv_opt] = fmincon(obj_fun, x_ga, [], [], [], [], lb, ub, @(x) constraints(x, ast), fm_opts);
    
    run_time = toc;
    
    %salvo nell'array delle run
    results_dv(k) = dv_opt;
    results_x(k, :) = x_opt;
    results_time(k) = run_time;
    
    if  dv_opt < best_dv_global
        best_dv_global = dv_opt;
        best_x_global = x_opt;
        best_history_ga = history_ga;
        best_history_fmincon = history_fmincon;
    end
end %fine run

% --- NUOVI CALCOLI STATISTICI ---
media_dv = mean(results_dv);
mediana_dv = median(results_dv);
varianza_dv = var(results_dv);
dev_std_dv = std(results_dv);
min_dv = min(results_dv);
max_dv = max(results_dv);
tempo_medio = mean(results_time);
tasso_successo = (length(results_dv) / N_runs) * 100;

% Calcolo dei percentili
prctile_25 = prctile(results_dv, 25);
prctile_75 = prctile(results_dv, 75);
iqr_dv = prctile_75 - prctile_25; % Interquartile Range

fprintf('\n===================================================\n');
fprintf('                REPORT STATISTICO                  \n');
fprintf('===================================================\n');
fprintf('Numero di run totali:    %d\n', N_runs);
fprintf('Tasso di successo:       %.1f%%\n', tasso_successo);
fprintf('Tempo computazionale m.: %.2f secondi/run\n', tempo_medio);
fprintf('---------------------------------------------------\n');
fprintf('Delta V Minimo (OTTIMO): %.4f km/s\n', min_dv);
fprintf('Delta V Massimo:         %.4f km/s\n', max_dv);
fprintf('Media:                   %.4f km/s\n', media_dv);
fprintf('Mediana:                 %.4f km/s\n', mediana_dv);
fprintf('25° Percentile:          %.4f km/s\n', prctile_25);
fprintf('75° Percentile:          %.4f km/s\n', prctile_75);
fprintf('Range Interquartile IQR: %.4f km/s\n', iqr_dv);
fprintf('Deviazione Standard:     %.4f km/s\n', dev_std_dv);
fprintf('===================================================\n');
fprintf('MIGLIOR SET DI VARIABILI TROVATO:\n');
fprintf('Theta_1 Partenza:        %.2f deg\n', rad2deg(best_x_global(1)));
fprintf('Theta_2 Arrivo:          %.2f deg\n', rad2deg(best_x_global(2)));
fprintf('Omega_T Trasferimento:   %.2f deg\n', rad2deg(best_x_global(3)));

% =========================================================================
% PLOTTING
% =========================================================================

figure('Name', 'Dispersione Statistica (Punti)');
hold on; grid on;

% Disegno tutti i punti (uno per ogni run)
scatter(1:N_runs, results_dv, 40, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeAlpha', 0.6, 'DisplayName', 'Singola Run');

% Evidenzio la run migliore in assoluto con una stella verde gigante
[~, best_run_idx] = min(results_dv);
plot(best_run_idx, min_dv, 'p', 'MarkerSize', 15, 'MarkerFaceColor', '#77AC30', 'MarkerEdgeColor', 'k', 'DisplayName', 'Ottimo Assoluto');

% Linee di riferimento orizzontali
yline(media_dv, '--r', 'LineWidth', 2, 'DisplayName', 'Media');
yline(prctile_75, ':k', 'LineWidth', 1.5, 'DisplayName', '75° Percentile');

title('Dispersione dei Costi per ogni singola Run (Monte Carlo)');
xlabel('Numero della Run');
ylabel('\DeltaV [km/s]');
legend('Location', 'northeast');

% --- FORZATURA ASSI (Zoom sui valori reali, ignora lo zero) ---
% Troviamo il range reale dei dati
min_reale = min(results_dv);
max_reale = max(results_dv);
delta_reale = max_reale - min_reale;

% Se l'algoritmo non ha mai sballato (range piccolo), aggiungiamo un margine fisso.
% Altrimenti usiamo un margine percentuale (5%) del range totale.
if delta_reale < 1e-5 % Caso molto compatto
    margine = 1e-6; 
else
    margine = 0.05 * delta_reale; % 5% di margine
end

% Impostiamo i limiti centrati sui dati
ylim([min_reale - margine, max_reale + margine]);
% -----------------------------------------------------------------

hold off;

% 2. Boxplot (Ottimo per individuare le run anomale)
figure('Name', 'Boxplot Analisi');
boxplot(results_dv, 'Labels', {'\DeltaV (Monte Carlo)'});
title('Boxplot: Dispersione e Outliers');
ylabel('\DeltaV [km/s]');
grid on;

% 3. Scatter Plot 3D dello Spazio delle Variabili (Mostra i minimi locali)
figure('Name', 'Mappa Variabili 3D');
scatter3(rad2deg(results_x(:,1)), rad2deg(results_x(:,2)), rad2deg(results_x(:,3)), ...
         60, results_dv, 'filled', 'MarkerEdgeColor', 'k');
colorbar;
title('Soluzioni Trovate nello Spazio di Ricerca');
xlabel('\theta_1 Partenza [deg]');
ylabel('\theta_2 Arrivo [deg]');
zlabel('\omega_T Trasferimento [deg]');
grid on; view(45, 30);

% 4. Funzione di Distribuzione Cumulativa (Probabilità di successo)
figure('Name', 'Probabilità Cumulativa');
dv_sorted = sort(results_dv);
prob = (1:N_runs) / N_runs;
stairs(dv_sorted, prob, 'LineWidth', 2.5, 'Color', '#77AC30');
grid on;
title('Funzione di Probabilità Cumulativa (ECDF)');
xlabel('\DeltaV [km/s]');
ylabel('Probabilità cumulativa');

% =========================================================================
% 5. PLOTTING DELLE CURVE DI CONVERGENZA (Run Migliore)
% =========================================================================
fprintf('\nGenerazione dei grafici di convergenza della run migliore...\n');

% --- FIGURE 5: Convergenza Algoritmo Genetico (GA) ---
if ~isempty(best_history_ga)
    figure('Name', 'Convergenza GA');
    
    pop_size = ga_opts.PopulationSize; 
    num_righe = size(best_history_ga, 1);
    num_gen = floor(num_righe / pop_size);
    
    best_dv_ga = zeros(num_gen, 1);
    mean_dv_ga = zeros(num_gen, 1);
    all_dv_ga = zeros(num_gen, pop_size); 
    
    for g = 1:num_gen
        idx_start = (g-1)*pop_size + 1;
        idx_end = g*pop_size;
        
        scores_gen = best_history_ga(idx_start:idx_end, 4);
        
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
    all_dv_ga(all_dv_ga >= 1e5) = NaN;
    
    h_altri = plot(1:num_gen, all_dv_ga, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    h_best = plot(1:num_gen, best_dv_ga, '-b', 'LineWidth', 2.5);
    h_mean = plot(1:num_gen, mean_dv_ga, '-c', 'LineWidth', 2);
    
    grid on;
    title('Fase 1: Convergenza Algoritmo Genetico');
    xlabel('Generazioni');
    ylabel('\DeltaV Totale [km/s]');
    legend([h_best, h_mean, h_altri(1)], ...
        {'Miglior Individuo', 'Media Popolazione', 'Altri Individui'}, 'Location', 'northeast');
end

% --- FIGURE 6: Convergenza fmincon ---
if ~isempty(best_history_fmincon)
    figure('Name', 'Convergenza fmincon');
    
    iters_fmincon = 1:size(best_history_fmincon, 1);
    dv_fmincon = best_history_fmincon(:, 4);
    
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