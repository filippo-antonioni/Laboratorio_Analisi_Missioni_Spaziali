function [interseca, th1_int, th2_int] = orbit_intersection(a1, e1, om1, a2, e2, om2)
    % Calcolo dei semilati retti
    p1 = a1 * (1 - e1^2);
    p2 = a2 * (1 - e2^2);
    
    % Calcolo dei coefficienti
    A = p2 * e1 * cos(om1) - p1 * e2 * cos(om2);
    B = p2 * e1 * sin(om1) - p1 * e2 * sin(om2);
    C = p1 - p2;
    
    % Condizione matematica di intersezione
    if A^2 + B^2 >= C^2
        interseca = true;
        
        phi = atan2(B, A);
        alpha = acos(C / sqrt(A^2 + B^2));
        
        % Angoli ASSOLUTI di intersezione (Argomenti di latitudine, u)
        u_int = [phi + alpha, phi - alpha];
        
        % 1. Trasformo nell'Anomalia Vera dell'Orbita 1 e 2
        th1_grezz = mod(u_int - om1, 2*pi); % Rispetto al pericentro 1
        th2_grezz = mod(u_int - om2, 2*pi); % Rispetto al pericentro 2
        
        % 2. Ordino i valori dell'Orbita 1 dal più piccolo al più grande
        [th1_int, idx_sort] = sort(th1_grezz);
        
        % 3. Riordino anche l'Orbita 2 usando gli stessi indici
        th2_int = th2_grezz(idx_sort);
        
    else
        interseca = false;
        th1_int = [];
        th2_int = [];
    end
end