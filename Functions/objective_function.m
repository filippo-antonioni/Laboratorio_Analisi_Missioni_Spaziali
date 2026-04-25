function dv = objective_function(x, ast)
    % Estraggo le variabili scelte dall'ottimizzatore in questa iterazione
    th1_i = x(1);   % Anomalia vera Terra (partenza)
    th2_f = x(2);   % Anomalia vera Asteroide (arrivo)
    om_T  = x(3);   % Argomento del pericentro del trasferimento (variabile libera)
    
    % --- DATI FISICI ---
    mu_sun = 1.32712440018e11; % Parametro gravitazionale Sole [km^3/s^2]
    R_sun = 696340;            % Raggio Sole [km]
    AU = 149597870.7;          % Unità Astronomica [km]
    
    % Parametri orbitali della Terra
    a_terra = 1.4946e8; 
    e_terra = 0.016; 
    i_terra = 9.1920e-5; 
    OM_terra = 2.7847; 
    om_terra = 5.2643;
    
    % =====================================================================
    % STEP 1: Calcolo Posizione e Velocità su orbita iniziale e finale
    % =====================================================================
    % NOTA: Assumo che tu abbia a disposizione la funzione par2car del Lab.
    % Se i vettori r1 e r2 sono allineati al 100%, il prodotto vettoriale fallisce.
    
    [r1, v1] = par2car(a_terra, e_terra, i_terra, OM_terra, om_terra, th1_i, mu_sun);
    [r2, v2] = par2car(ast.a, ast.e, ast.i, ast.OM, ast.om, th2_f, mu_sun);
    
    % =====================================================================
    % STEP 2: Definizione del piano di trasferimento
    % =====================================================================
    h_vect = cross(r1, r2);
    norm_h = norm(h_vect);
    
    if norm_h < 1e-6
        dv = 1e6; return; 
    end
    
    h_T = h_vect / norm_h; % Versore normale al piano
    
    i_T = acos(h_T(3)); % Inclinazione
    
    % Calcolo della Linea dei Nodi e RAAN
    K = [0; 0; 1];
    N_vect = cross(K, h_T);
    N_T = N_vect / norm(N_vect);
    
    if N_T(2) >= 0
        OM_T = acos(N_T(1));
    else
        OM_T = 2*pi - acos(N_T(1));
    end
    
    % =====================================================================
    % STEP 3: Sistema Perifocale e Anomalie di Trasferimento
    % =====================================================================
    % Matrici di rotazione per passare da Inerziale a Perifocale
    R_OM = [cos(OM_T) sin(OM_T) 0; -sin(OM_T) cos(OM_T) 0; 0 0 1];
    R_i  = [1 0 0; 0 cos(i_T) sin(i_T); 0 -sin(i_T) cos(i_T)];
    R_om = [cos(om_T) sin(om_T) 0; -sin(om_T) cos(om_T) 0; 0 0 1];
    
    T_matrix = R_om * R_i * R_OM; 
    
    % Posizioni nel sistema perifocale
    r1_PF = T_matrix * r1;
    r2_PF = T_matrix * r2;
    
    % Anomalie vere sull'orbita di trasferimento
    th1_T = atan2(r1_PF(2), r1_PF(1));
    th2_T = atan2(r2_PF(2), r2_PF(1));
    
    % =====================================================================
    % STEP 4: Forma dell'orbita di Trasferimento (a_T, e_T)
    % =====================================================================
    r1_mag = norm(r1);
    r2_mag = norm(r2);
    
    % Ricavo l'eccentricità eguagliando l'equazione della conica
    num_e = r2_mag - r1_mag;
    den_e = r1_mag * cos(th1_T) - r2_mag * cos(th2_T);
    
    if abs(den_e) < 1e-6
        dv = 1e6; return;
    end
    
    e_T = num_e / den_e;
    
    % Se l'eccentricità non è fisicamente valida (negativa o non chiusa), scarta!
    if e_T < 0 || e_T >= 0.999
        dv = 1e6; return;
    end
    
    % Calcolo il semiasse maggiore
    a_T = (r1_mag * (1 + e_T * cos(th1_T))) / (1 - e_T^2);
    
    % Controllo di sicurezza: non schiantarsi sul Sole
    rp_T = a_T * (1 - e_T);
    if rp_T < (R_sun + 50000000) % Aggiungo 50.000.000 km di margine
        dv = 1e6; return;
    end
    
    % =====================================================================
    % STEP 5: Calcolo Velocità e Delta V [Slide 15]
    % =====================================================================
    % Ora che ho [a_T, e_T, i_T, OM_T, om_T], uso di nuovo par2car per trovare 
    % le velocità ai punti 1 e 2 sull'orbita di trasferimento [cite: 212-214]
    
    [~, v1_T] = par2car(a_T, e_T, i_T, OM_T, om_T, th1_T, mu_sun);
    [~, v2_T] = par2car(a_T, e_T, i_T, OM_T, om_T, th2_T, mu_sun);
    
    % Calcolo il Delta V totale vettoriale [cite: 217]
    dV1 = norm(v1_T - v1);
    dV2 = norm(v2 - v2_T);
    
    dv = dV1 + dV2; % Questo è il vero costo in km/s!
end