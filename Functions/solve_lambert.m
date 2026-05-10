function [v1, v2] = solve_lambert(r1, r2, dt, mu, prograde)
    % SOLVE_LAMBERT Risolutore del Problema di Lambert (Variabili Universali)
    %
    % Input:
    %   r1       - vettore posizione iniziale (3x1) [km]
    %   r2       - vettore posizione finale (3x1) [km]
    %   dt       - tempo di volo [s]
    %   mu       - parametro gravitazionale del corpo centrale [km^3/s^2]
    %   prograde - 1 per orbita prograda (default), 0 per retrograda
    %
    % Output:
    %   v1       - vettore velocità iniziale (3x1) [km/s]
    %   v2       - vettore velocità finale (3x1) [km/s]

    if nargin < 5
        prograde = 1; % Di default assume orbita prograda
    end

    % Assicuriamoci che r1 e r2 siano vettori colonna
    r1 = r1(:);
    r2 = r2(:);

    r1_norm = norm(r1);
    r2_norm = norm(r2);

    % 1. Calcolo del delta theta (angolo di trasferimento)
    cross12 = cross(r1, r2);
    theta = acos(dot(r1, r2) / (r1_norm * r2_norm));
    
    % Controllo per decidere il verso di percorrenza
    if prograde == 1
        if cross12(3) <= 0
            theta = 2*pi - theta;
        end
    else
        if cross12(3) >= 0
            theta = 2*pi - theta;
        end
    end

    % 2. Calcolo della costante geometrica A
    A = sin(theta) * sqrt(r1_norm * r2_norm / (1 - cos(theta)));

    % 3. Newton-Raphson per trovare la variabile universale z
    z = 0; % Ipotesi iniziale (orbita parabolica)
    tol = 1e-6;
    n_max = 500;
    ratio = 1;
    
    for n = 1:n_max
        if abs(ratio) < tol
            break;
        end
        
        [C, S] = stumpff(z);
        y = r1_norm + r2_norm + A * (z*S - 1) / sqrt(C);
        
        % Controllo di sicurezza: se y diventa negativo, forza un valore piccolo 
        % per evitare che la radice quadrata restituisca numeri complessi
        if A > 0 && y < 0
            y = 1e-5; 
        end
        
        chi = sqrt(y / C);
        t_calc = (chi^3 * S + A * sqrt(y)) / sqrt(mu);
        
        % Derivata dt/dz calcolata numericamente (più stabile di quella analitica)
        dz = 1e-5;
        [C_d, S_d] = stumpff(z + dz);
        y_d = r1_norm + r2_norm + A * ((z+dz)*S_d - 1) / sqrt(C_d);
        if A > 0 && y_d < 0; y_d = 1e-5; end
        chi_d = sqrt(y_d / C_d);
        t_calc_d = (chi_d^3 * S_d + A * sqrt(y_d)) / sqrt(mu);
        
        dt_dz = (t_calc_d - t_calc) / dz;
        
        % Aggiornamento di z
        ratio = (t_calc - dt) / dt_dz;
        z = z - ratio;
    end
    
    if n >= n_max
        warning('Lambert: Newton-Raphson non è arrivato a convergenza ottimale.');
    end

    % 4. Calcolo dei Coefficienti di Lagrange (f, g, f_dot, g_dot)
    [C, S] = stumpff(z);
    y = r1_norm + r2_norm + A * (z*S - 1) / sqrt(C);
    
    f = 1 - y / r1_norm;
    g = A * sqrt(y / mu);
    g_dot = 1 - y / r2_norm;

    % 5. Calcolo dei vettori velocità in 3D
    v1 = (1/g) * (r2 - f * r1);
    v2 = (1/g) * (g_dot * r2 - r1);
end

% --- FUNZIONI AUSILIARIE ---
function [C, S] = stumpff(z)
    % Calcola le funzioni di Stumpff C(z) e S(z)
    if z > 0
        sqrt_z = sqrt(z);
        S = (sqrt_z - sin(sqrt_z)) / (sqrt_z^3);
        C = (1 - cos(sqrt_z)) / z;
    elseif z < 0
        sqrt_z = sqrt(-z);
        S = (sinh(sqrt_z) - sqrt_z) / (sqrt_z^3);
        C = (cosh(sqrt_z) - 1) / (-z);
    else
        S = 1/6;
        C = 1/2;
    end
end