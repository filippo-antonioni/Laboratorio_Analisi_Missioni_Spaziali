clear
close all
clc

% --- INIZIALIZZAZIONE STORICO ---
global history_pso history_fmincon;
history_pso = [];    
history_fmincon = [];

mu = 398600;

% --- Parametri orbita eliocentrica iniziale (Terra) ---
a_i=1.4946*1e8; 
e_i=0.016; 
i_i=9.1920*1e-5; 
OM_i=2.7847; 
om_i=5.2643;

% --- Dati dell'asteroide 363505 (2003 UC20) ---
AU_to_km = 149597870.7; 
ast.a = 0.781241 * AU_to_km; 
ast.e = 0.336932; 
ast.i = deg2rad(3.78); 
ast.OM = deg2rad(187.92); 
ast.om = deg2rad(60.16);

% 2. CONFIGURAZIONE OTTIMIZZAZIONE
lb = [0, 0, 0]; ub = [2*pi, 2*pi, 2*pi];
obj_fun = @(x) objective_function(x, ast);

% --- FASE 1: RICERCA GLOBALE CON PSO ---
fprintf('Fase 1: Avvio Particle Swarm Optimization (PSO)...\n');
pso_opts = optimoptions('particleswarm', 'SwarmSize', 200, 'MaxIterations', 100, ...
    'Display', 'iter', 'OutputFcn', @outfun_pso);

[x_pso, fval_pso] = particleswarm(obj_fun, 3, lb, ub, pso_opts);

% --- FASE 2: RIFINITURA LOCALE CON FMINCON --- 
fprintf('\nFase 2: Rifinitura con fmincon (Partendo dal risultato PSO)...\n');
fm_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'TolFun', 1e-8, 'OutputFcn', @outfun_fmincon_pso);

[x_opt, dv_opt] = fmincon(obj_fun, x_pso, [], [], [], [], lb, ub, @(x) constraints(x, ast), fm_opts);

% 3. RISULTATI E PLOTTING
fprintf('\n--- RISULTATI OTTIMI ---\nDelta V Totale: %.4f km/s\n', dv_opt);

% =========================================================================
% 4. PLOTTING DELLE CURVE DI CONVERGENZA
% =========================================================================

% --- FIGURE 1: Convergenza PSO ---
if ~isempty(history_pso)
    figure('Name', 'Convergenza PSO', 'Color', 'w', 'Position', [100, 100, 700, 500]);
    hold on; grid on;
    
    iterazioni = history_pso(:,1);
    best_val_pso = history_pso(:,2);
    swarm_vals = history_pso(:, 3:end);
    
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
    
    title('Fase 1: Convergenza Particle Swarm (PSO)');
    xlabel('Iterazioni'); ylabel('\DeltaV Totale [km/s]');
    
    % Legenda aggiornata
    legend([h_best, h_mean, h_altri(1)], ...
        {'Miglior Particella (Global Best)', 'Media dello Sciame', 'Altre Particelle dello Sciame'}, ...
        'Location', 'northeast');
end

% --- FIGURE 2: Convergenza fmincon ---
if ~isempty(history_fmincon)
    figure('Name', 'Convergenza fmincon', 'Color', 'w', 'Position', [150, 150, 700, 500]);
    plot(history_fmincon(:,1), history_fmincon(:,2), '-o', 'Color', '#D95319', 'LineWidth', 1.5, 'MarkerFaceColor', '#EDB120');
    grid on;
    title('Fase 2: Rifinitura fmincon');
    xlabel('Iterazioni'); ylabel('\DeltaV Totale [km/s]');
end