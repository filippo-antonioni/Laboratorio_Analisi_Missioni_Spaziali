function [c, ceq] = constraints(x, ast)
    % Estraggo le variabili scelte dall'ottimizzatore
    th1_i = x(1);   
    th2_f = x(2);   
    om_T  = x(3);   
    
    % --- DATI FISICI ---
    mu_sun = 1.32712440018e11; 
    R_sun = 696340; 
    
    % Parametri orbitali della Terra
    a_terra = 1.4946e8; 
    e_terra = 0.016; 
    i_terra = 9.1920e-5; 
    OM_terra = 2.7847; 
    om_terra = 5.2643;
    
    % --- RICALCOLO eT e aT PER IL TENTATIVO CORRENTE x ---
    [r1, ~] = par2car(a_terra, e_terra, i_terra, OM_terra, om_terra, th1_i, mu_sun);
    [r2, ~] = par2car(ast.a, ast.e, ast.i, ast.OM, ast.om, th2_f, mu_sun);
    
    h_vect = cross(r1, r2);
    norm_h = norm(h_vect);
    
    if norm_h < 1e-6
        c = [1; 1; 1]; ceq = []; return; % Forza violazione per vettori allineati
    end
    
    h_T = h_vect / norm_h; 
    i_T = acos(h_T(3)); 
    
    K = [0; 0; 1];
    N_vect = cross(K, h_T);
    norm_N = norm(N_vect);
    
    if norm_N < 1e-6
        OM_T = 0;
    else
        N_T = N_vect / norm_N;
        if N_T(2) >= 0
            OM_T = acos(N_T(1));
        else
            OM_T = 2*pi - acos(N_T(1));
        end
    end
    
    R_OM = [cos(OM_T) sin(OM_T) 0; -sin(OM_T) cos(OM_T) 0; 0 0 1];
    R_i  = [1 0 0; 0 cos(i_T) sin(i_T); 0 -sin(i_T) cos(i_T)];
    R_om = [cos(om_T) sin(om_T) 0; -sin(om_T) cos(om_T) 0; 0 0 1];
    T_matrix = R_om * R_i * R_OM; 
    
    r1_PF = T_matrix * r1;
    r2_PF = T_matrix * r2;
    
    th1_T = atan2(r1_PF(2), r1_PF(1));
    th2_T = atan2(r2_PF(2), r2_PF(1));
    
    r1_mag = norm(r1);
    r2_mag = norm(r2);
    
    num_e = r2_mag - r1_mag;
    den_e = r1_mag * cos(th1_T) - r2_mag * cos(th2_T);
    
    if abs(den_e) < 1e-6
        c = [1; 1; 1]; ceq = []; return; % Forza violazione per orbita degenere
    end
    
    eT = num_e / den_e;
    aT = (r1_mag * (1 + eT * cos(th1_T))) / (1 - eT^2);
    rpT = aT * (1 - eT);
    
    % --- VINCOLI ---
    c(1) = eT - 0.999;               % Deve essere ellittica chiusa
    c(2) = -eT;                      % Eccentricità >= 0
    c(3) = (R_sun + 50000000) - rpT; % Non colpire il Sole
    
    ceq = []; % Nessun vincolo di uguaglianza
end