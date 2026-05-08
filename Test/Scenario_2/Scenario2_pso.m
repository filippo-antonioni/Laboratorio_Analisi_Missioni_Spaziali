clear
close all
clc

% --- INIZIALIZZAZIONE STORICO ---
global history_pso history_fmincon;
history_pso = [];    
history_fmincon = [];

best_dv_global = inf;
best_x_global = [];
best_history_pso = [];
best_history_fmincon = [];

N_runs = 200;
results_x = zeros(N_runs, 3);
results_dv = zeros(N_runs, 1);
results_time = zeros(N_runs, 1);
mu = 398600;

% --- Parametri orbita eliocentrica iniziale (Terra) ---
a_i = 1.4946*1e8; 
e_i = 0.016; 
i_i = 9.1920*1e-5; 
OM_i = 2.7847; 
om_i = 5.2643;

% --- Dati dell'asteroide 363505 (2003 UC20) ---
AU_to_km = 149597870.7; 
ast.a = 0.781241 * AU_to_km; 
ast.e = 0.336932; 
ast.i = deg2rad(3.78); 
ast.OM = deg2rad(187.92); 
ast.om = deg2rad(60.16);

% 2. CONFIGURAZIONE OTTIMIZZAZIONE
% Variabili x = [theta1_partenza, theta2_arrivo, omega_trasferimento]
lb = [0, 0, 0]; 
ub = [2*pi, 2*pi, 2*pi];
obj_fun = @(x) objective_function(x, ast);

fprintf('--- AVVIO ANALISI STATISTICA CON PSO (%d RUNS) ---\n', N_runs);

% --- CICLO DI OTTIMIZZAZIONE MONTE CARLO ---
for k = 1:N_runs
    
    
    % Azzero le variabili globali per la singola run
    history_pso = [];
    history_fmincon = [];
    
    tic % Inizio a contare il tempo
    
    % --- FASE 1: RICERCA GLOBALE CON PSO ---
    pso_opts = optimoptions('particleswarm', 'SwarmSize', 200, 'MaxIterations', 100, ...
        'Display', 'off', 'OutputFcn', @outfun_pso);
    [x_pso, fval_pso] = particleswarm(obj_fun, 3, lb, ub, pso_opts);
    
    % --- FASE 2: RIFINITURA LOCALE CON FMINCON --- 
    fm_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'TolFun', 1e-8, 'OutputFcn', @outfun_fmincon_pso);
    [x_opt, dv_opt] = fmincon(obj_fun, x_pso, [], [], [], [], lb, ub, @(x) constraints(x, ast), fm_opts);
    
    run_time = toc; % Fine cronometro
    
    % Salvo nell'array delle run
    results_dv(k) = dv_opt;
    results_x(k, :) = x_opt;
    results_time(k) = run_time;
    
    
    
    % Aggiorno il record assoluto e salvo la storia per i plot finali
    if dv_opt < best_dv_global
        best_dv_global = dv_opt;
        best_x_global = x_opt;
        best_history_pso = history_pso;
        best_history_fmincon = history_fmincon;
    end
end % fine run

% --- CALCOLI STATISTICI ---
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
% PLOTTING STATISTICO
% =========================================================================

% 1. Scatter Plot della dispersione dei Delta V (Punti Singoli centrati)
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

title('Dispersione dei Costi per ogni singola Run (Monte Carlo - PSO)');
xlabel('Numero della Run');
ylabel('\DeltaV [km/s]');
legend('Location', 'northeast');

% --- FORZATURA ASSI (Zoom sui valori reali, ignora lo zero) ---
min_reale = min(results_dv);
max_reale = max(results_dv);
delta_reale = max_reale - min_reale;

if delta_reale < 1e-5 % Caso molto compatto
    margine = 1e-6; 
else
    margine = 0.05 * delta_reale; % 5% di margine
end
ylim([min_reale - margine, max_reale + margine]);
hold off;

% 2. Boxplot (Ottimo per individuare le run anomale)
figure('Name', 'Boxplot Analisi');
boxplot(results_dv, 'Labels', {'\DeltaV (Monte Carlo)'});
title('Boxplot: Dispersione e Outliers (PSO)');
ylabel('\DeltaV [km/s]');
grid on;

% 3. Scatter Plot 3D dello Spazio delle Variabili (Mostra i minimi locali)
figure('Name', 'Mappa Variabili 3D');
scatter3(rad2deg(results_x(:,1)), rad2deg(results_x(:,2)), rad2deg(results_x(:,3)), ...
         60, results_dv, 'filled', 'MarkerEdgeColor', 'k');
colorbar;
title('Soluzioni Trovate nello Spazio di Ricerca (PSO)');
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

% --- FIGURE 5: Convergenza PSO ---
if ~isempty(best_history_pso)
    figure('Name', 'Convergenza PSO');
    hold on; grid on;
    
    iterazioni = best_history_pso(:,1);
    best_val_pso = best_history_pso(:,2);
    swarm_vals = best_history_pso(:, 3:end);
    
    % Pulizia valori per il plot (evitiamo che i 1e6 schiaccino il grafico)
    swarm_vals(swarm_vals > 100) = NaN; 
    
    % Calcolo della media dello sciame ad ogni iterazione
    mean_val_pso = mean(swarm_vals, 2, 'omitnan');
    
    % 1. Plot di tutte le particelle in grigio chiaro
    h_altri = plot(iterazioni, swarm_vals, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    
    % 2. Plot della Media in ciano (linea continua)
    h_mean = plot(iterazioni, mean_val_pso, '-c', 'LineWidth', 2);
    
    % 3. Plot del Global Best in blu (linea continua spessa)
    h_best = plot(iterazioni, best_val_pso, '-b', 'LineWidth', 2.5);
    
    title('Fase 1: Convergenza Particle Swarm (Run Migliore)');
    xlabel('Iterazioni'); ylabel('\DeltaV Totale [km/s]');
    
    % Legenda aggiornata
    legend([h_best, h_mean, h_altri(1)], ...
        {'Miglior Particella (Global Best)', 'Media dello Sciame', 'Altre Particelle dello Sciame'}, ...
        'Location', 'northeast');
end

% --- FIGURE 6: Convergenza fmincon ---
if ~isempty(best_history_fmincon)
    figure('Name', 'Convergenza fmincon');
    plot(best_history_fmincon(:,1), best_history_fmincon(:,2), '-o', 'Color', '#D95319', 'LineWidth', 1.5, 'MarkerFaceColor', '#EDB120');
    grid on;
    title('Fase 2: Rifinitura fmincon (Run Migliore)');
    xlabel('Iterazioni'); ylabel('\DeltaV Totale [km/s]');
end