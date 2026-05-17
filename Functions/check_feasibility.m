function [is_feasible, reason] = check_feasibility(r_p, r_a, R_ast, margine)
    % CHECK_FEASIBILITY Verifica la validità geometrica di un'orbita attorno a un corpo.
    %
    % INPUT:
    % r_p      - Raggio di pericentro calcolato [km]
    % r_a      - Raggio di apocentro calcolato [km]
    % R_ast    - Raggio fisico dell'asteroide/pianeta [km]
    % R_SOI    - Raggio della Sfera di Influenza [km]
    % margine  - (Opzionale) Quota di sicurezza minima dalla superficie [km]
    %
    % OUTPUT:
    % is_feasible - true se l'orbita è sicura, false se non lo è
    % reason      - Stringa che spiega il motivo del fallimento
    
    % Se non passiamo il margine alla funzione, impostiamo un default (es. 100 metri)
    if nargin < 4
        margine = 0.1; 
    end
    
    is_feasible = true;
    reason = 'OK';
    
    % 1. Controllo Impatto: Il pericentro è sotto o sulla superficie?
    if r_p <= R_ast
        is_feasible = false;
        reason = 'IMPATTO: Il pericentro è sotto la superficie dell''asteroide.';
        return; % Esce immediatamente dalla funzione
    end
    
    % 2. Controllo Quota di Sicurezza: Passa troppo vicino?
    if r_p < (R_ast + margine)
        is_feasible = false;
        reason = 'TROPPO VICINO: Orbita sotto il margine di sicurezza impostato.';
        return;
    end
    
    % 3. Controllo Coerenza Geometrica: L'apocentro è minore del pericentro?
    if r_a < r_p
        is_feasible = false;
        reason = 'NON FISICO: L''apocentro calcolato è minore del pericentro.';
        return;
    end

end