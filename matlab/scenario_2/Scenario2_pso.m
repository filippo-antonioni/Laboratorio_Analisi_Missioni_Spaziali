clear
close all
clc

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

% Parametri orbita eliocentrica della Terra
a_i = 1.4946*1e8; 
e_i = 0.016; 
i_i = 9.1920*1e-5; 
OM_i = 2.7847; 
om_i = 5.2643;
mu_sun = 1.32712440018e11;

% Dati dell'asteroide 363505 (2003 UC20)
AU_to_km = 149597870.7; 
ast.a = 0.781241 * AU_to_km; 
ast.e = 0.336932; 
ast.i = deg2rad(3.78); 
ast.OM = deg2rad(187.92); 
ast.om = deg2rad(60.16);

% Variabili x = [theta1_partenza, theta2_arrivo, omega_trasferimento]
lb = [0, 0, 0]; 
ub = [2*pi, 2*pi, 2*pi];
obj_fun = @(x) objective_function(x, ast);

% RICERCA PARAMETRI PSO
toll = 1e-4; 
max_iter = 20; 
iter = 0;
swarm_size = 0;  
max_iterations = 0; 
N_runs_pop = 80; 
err = inf; 
storia_swarm = [];
storia_err = [];

while iter < max_iter && err > toll 
    iter = iter + 1;
    swarm_size = swarm_size + 10; % Aumento combinato
    max_iterations = max_iterations + 10;
    
    results_dv_temp = zeros(N_runs_pop, 1);
    
    for k = 1:N_runs_pop
        % Fase 1: PSO
        pso_opts = optimoptions('particleswarm', 'SwarmSize', swarm_size, 'MaxIterations', max_iterations, 'Display', 'off');
        [x_pso, ~] = particleswarm(obj_fun, 3, lb, ub, pso_opts);
        
        % Fase 2: fmincon
        fm_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
        [~, dv_opt] = fmincon(obj_fun, x_pso, [], [], [], [], lb, ub, @(x) constraints(x, ast), fm_opts);
        
        results_dv_temp(k) = dv_opt;
    end 
    
    err = std(results_dv_temp); 
  
    % Salvo i dati
    storia_swarm = [storia_swarm, swarm_size];
    storia_err = [storia_err, err];
    
   
end

if err <= toll
    fprintf('\nConvergenza raggiunta: SwarmSize=%d, MaxIter=%d\n', swarm_size, max_iterations);
else
    fprintf('\nConvergenza raggiunta: SwarmSize=%d, MaxIter=%d\n', swarm_size, max_iterations);
end

% PLOT ANDAMENTO ERRORE PSO
figure('Name', 'Andamento Errore PSO');
plot(storia_swarm, storia_err, '-o', 'LineWidth', 2.5, 'MarkerFaceColor', '#D95319', 'MarkerSize', 8, 'Color', '#D95319');
hold on;
yline(toll, '--r', 'LineWidth', 2, 'DisplayName', 'Tolleranza Target');
grid on;
if max(storia_err) > 10
    ylim([0, 10]); 
end
title('Crollo della Dispersione all''aumentare di SwarmSize/MaxIterations');
xlabel('SwarmSize e MaxIterations (Valore)');
ylabel('Deviazione Standard (Errore) [km/s]');
legend('Errore misurato', 'Tolleranza Target', 'Location', 'northeast');
hold off;

% Fissiamo i parametri trovati per il ciclo finale
opt_swarm_size = swarm_size;
opt_max_iterations = max_iterations;

fprintf('\n--- AVVIO ANALISI STATISTICA CON PSO (%d RUNS) ---\n', N_runs);

for k = 1:N_runs
    
    history_pso = [];
    history_fmincon = [];
    
    tic % Inizio a contare il tempo
    
    % RICERCA CON PSO
    pso_opts = optimoptions('particleswarm', 'SwarmSize', opt_swarm_size, 'MaxIterations', opt_max_iterations, ...
        'Display', 'off', 'OutputFcn', @outfun_pso);
    [x_pso, fval_pso] = particleswarm(obj_fun, 3, lb, ub, pso_opts);
    
    % FMINCON 
    fm_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'TolFun', 1e-8, 'OutputFcn', @outfun_fmincon_pso);
    [x_opt, dv_opt] = fmincon(obj_fun, x_pso, [], [], [], [], lb, ub, @(x) constraints(x, ast), fm_opts);
    
    run_time = toc;
    
    results_dv(k) = dv_opt;
    results_x(k, :) = x_opt;
    results_time(k) = run_time;
    
    if dv_opt < best_dv_global
        best_dv_global = dv_opt;
        best_x_global = x_opt;
        best_history_pso = history_pso;
        best_history_fmincon = history_fmincon;
    end
end % fine run

media_dv = mean(results_dv);
mediana_dv = median(results_dv);
varianza_dv = var(results_dv);
dev_std_dv = std(results_dv);
min_dv = min(results_dv);
max_dv = max(results_dv);
tempo_medio = mean(results_time);
tasso_successo = (length(results_dv) / N_runs) * 100;

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
fprintf('Deviazione Standard:     %.4f km/s\n', dev_std_dv);
fprintf('===================================================\n');
fprintf('MIGLIOR SET DI VARIABILI TROVATO:\n');
fprintf('Theta_1 Partenza:        %.2f deg\n', rad2deg(best_x_global(1)));
fprintf('Theta_2 Arrivo:          %.2f deg\n', rad2deg(best_x_global(2)));
fprintf('Omega_T Trasferimento:   %.2f deg\n', rad2deg(best_x_global(3)));

opt_th1 = best_x_global(1);
opt_th2 = best_x_global(2);
opt_omT = best_x_global(3);

% 1. Ricalcolo Raggi Vettore alla partenza e all'arrivo
[r1_opt, ~] = par2car(a_i, e_i, i_i, om_i, om_i, opt_th1, mu_sun);
[r2_opt, ~] = par2car(ast.a, ast.e, ast.i, ast.OM, ast.om, opt_th2, mu_sun);

% 2. Inclinazione (i_trasf_opt) e RAAN (OM_transf_opt)
h_vec_opt = cross(r1_opt, r2_opt);
h_vers_opt = h_vec_opt / norm(h_vec_opt);
i_trasf_opt = acos(h_vers_opt(3));

N_vec_opt = cross([0; 0; 1]', h_vers_opt')'; 
if norm(N_vec_opt) < 1e-6
    OM_transf_opt = 0;
else
    N_vers_opt = N_vec_opt / norm(N_vec_opt);
    if N_vers_opt(2) >= 0
        OM_transf_opt = acos(N_vers_opt(1));
    else
        OM_transf_opt = 2*pi - acos(N_vers_opt(1));
    end
end

R_OM_tmp = [ cos(OM_transf_opt),  sin(OM_transf_opt), 0;
            -sin(OM_transf_opt),  cos(OM_transf_opt), 0;
                   0,               0,        1];
R_i_tmp  = [1,        0,               0;
            0,  cos(i_trasf_opt),  sin(i_trasf_opt);
            0, -sin(i_trasf_opt),  cos(i_trasf_opt)];
R_om_tmp = [ cos(opt_omT),  sin(opt_omT), 0;
            -sin(opt_omT),  cos(opt_omT), 0;
                   0,         0,  1];

T_Elio_PF_tmp = R_om_tmp * R_i_tmp * R_OM_tmp;
r1_PF_tmp = T_Elio_PF_tmp * r1_opt;
r2_PF_tmp = T_Elio_PF_tmp * r2_opt;

th1_T_opt = atan2(r1_PF_tmp(2), r1_PF_tmp(1));
th2_T_opt = atan2(r2_PF_tmp(2), r2_PF_tmp(1));

norm_r1 = norm(r1_opt);
norm_r2 = norm(r2_opt);
den_e_tmp = norm_r1 * cos(th1_T_opt) - norm_r2 * cos(th2_T_opt);
e_trasf_opt = (norm_r2 - norm_r1) / den_e_tmp;
p_T_opt = norm_r1 * (1 + e_trasf_opt * cos(th1_T_opt));
a_T_opt = p_T_opt / (1 - e_trasf_opt^2);

fprintf('---------------------------------------------------\n');
fprintf('PARAMETRI ORBITA DI TRASFERIMENTO OTTIMA:\n');
fprintf('Inclinazione (i):        %.4f rad (%.2f deg)\n', i_trasf_opt, rad2deg(i_trasf_opt));
fprintf('RAAN (OM):               %.4f rad (%.2f deg)\n', OM_transf_opt, rad2deg(OM_transf_opt));
fprintf('Eccentricita (e):        %.4f\n', e_trasf_opt);
fprintf('===================================================\n');

% =========================================================================
% PLOT
% =========================================================================
% Funzione di Distribuzione Cumulativa (Probabilità di successo)
figure('Name', 'Probabilità Cumulativa');
dv_sorted = sort(results_dv);
prob = (1:N_runs) / N_runs;
stairs(dv_sorted, prob, 'LineWidth', 2.5, 'Color', '#77AC30');
grid on;
title('Funzione di Probabilità Cumulativa (ECDF)');
xlabel('\DeltaV [km/s]');
ylabel('Probabilità cumulativa');

% CURVE DI CONVERGENZA (Run Migliore)
fprintf('\nGenerazione dei grafici di convergenza della run migliore...\n');

% Convergenza PSO
if ~isempty(best_history_pso)
    figure('Name', 'Convergenza PSO');
    hold on; 
    grid on;
    
    iterazioni = best_history_pso(:,1);
    best_val_pso = best_history_pso(:,2);
    swarm_vals = best_history_pso(:, 3:end);
    
    swarm_vals(swarm_vals > 100) = NaN; 
    mean_val_pso = mean(swarm_vals, 2, 'omitnan');
    
    h_altri = plot(iterazioni, swarm_vals, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    h_mean = plot(iterazioni, mean_val_pso, '-c', 'LineWidth', 2);
    h_best = plot(iterazioni, best_val_pso, '-b', 'LineWidth', 2.5);
    
    title('Fase 1: Convergenza Particle Swarm (Run Migliore)');
    xlabel('Iterazioni'); ylabel('\DeltaV Totale [km/s]');
    
    legend([h_best, h_mean, h_altri(1)], ...
        {'Miglior Particella (Global Best)', 'Media dello Sciame', 'Altre Particelle'}, ...
        'Location', 'northeast');
end

% Convergenza fmincon
if ~isempty(best_history_fmincon)
    figure('Name', 'Convergenza fmincon');
    plot(best_history_fmincon(:,1), best_history_fmincon(:,2), '-o', 'Color', '#D95319', 'LineWidth', 1.5, 'MarkerFaceColor', '#EDB120');
    grid on;
    title('Fase 2: Rifinitura fmincon (Run Migliore)');
    xlabel('Iterazioni'); ylabel('\DeltaV Totale [km/s]');
end

% TOF
TOF_sec = TOF(a_T_opt, e_trasf_opt, th1_T_opt, th2_T_opt, mu_sun);
TOF_giorni = TOF_sec / (24 * 3600);