function [c, ceq] = constraints(x, ast)
    % Ricalcolo eT e aT per il tentativo corrente x
    % ... (stessa logica della funzione obiettivo) ...
    
    eT = 0.5; % Esempio
    rpT = 1.2e8; % Esempio raggio pericentro
    R_sun = 696340; % [km]
    
    % Vincoli di disuguaglianza (c <= 0)
    c(1) = eT - 0.999;      % Deve essere ellittica [cite: 159]
    c(2) = (R_sun + 5e7) - rpT; % Non colpire il Sole [cite: 161]
    
    ceq = []; % Nessun vincolo di uguaglianza
end