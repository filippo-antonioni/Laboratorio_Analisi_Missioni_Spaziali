function plotTransfer(state1, type1, state2, type2, maneuver, mu, varargin)
% plotTransfer Plot 3D di un trasferimento orbitale
%
% plotTransfer(state1, type1, state2, type2, maneuver, mu, varargin)
%
% -------------------------------------------------------------------------
% Input arguments:
% state1      [1x6] o [6x1]   Vettore di stato dell'orbita 1 
% type1       [char]          'car' (Cartesiano [x,y,z,vx,vy,vz]) o 
%                             'par' (Kepleriano [a, e, i, OM, om, th])
% state2      [1x6] o [6x1]   Vettore di stato dell'orbita 2 
% type2       [char]          'car' o 'par'
% maneuver    [char]          Tipo di manovra. Scegliere tra:
%                             'pa', 'ap', 'pp', 'aa' (Bitangenti)
%                             'bielliptic' (Richiede varargin{1} = ra_t)
%                             'plane_change'
%                             'pericenter_change'
% mu          [1x1]           Parametro gravitazionale planetario [km^3/s^2]
% varargin    [cell]          Parametri aggiuntivi (es. ra_t per biellittico)
% -------------------------------------------------------------------------

    % --- 1. Uniformazione dei Dati (Conversione in Kepleriani) ---
    if strcmpi(type1, 'car')
        [a1, e1, i1, OM1, om1, ~] = car2par(state1(1:3), state1(4:6), mu);
    else
        a1 = state1(1); e1 = state1(2); i1 = state1(3); OM1 = state1(4); om1 = state1(5);
    end
    
    if strcmpi(type2, 'car')
        [a2, e2, i2, OM2, om2, ~] = car2par(state2(1:3), state2(4:6), mu);
    else
        a2 = state2(1); e2 = state2(2); i2 = state2(3); OM2 = state2(4); om2 = state2(5);
    end
    
    % --- 2. Setup della Figura ---
    %figure('Name', 'Orbital Transfer 3D', 'NumberTitle', 'off');
    hold on; grid on; axis equal;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    
    dth = 0.05; % Step per l'anomalia vera
    
    % Funzione anonima per ottenere i punti della traiettoria 
    % (sfrutta la funzione annidata a fine script)
    get_traj = @(a, e, i, OM, om, th_start, th_end) calc_traj(a, e, i, OM, om, th_start, th_end, dth, mu);
    
    % --- 3. Plot delle Orbite Iniziale e Finale (Intere) ---
    traj1 = get_traj(a1, e1, i1, OM1, om1, 0, 2*pi);
    traj2 = get_traj(a2, e2, i2, OM2, om2, 0, 2*pi);
    
    plot3(traj1(1,:), traj1(2,:), traj1(3,:), 'b', 'LineWidth', 1.5, 'DisplayName', 'Orbita Iniziale');
    plot3(traj2(1,:), traj2(2,:), traj2(3,:), 'g', 'LineWidth', 1.5, 'DisplayName', 'Orbita Finale');
    
    % --- 4. Calcolo e Plot dell'Orbita di Trasferimento ---
    switch lower(maneuver)
        case {'pa', 'ap', 'pp', 'aa'}
            % TRASFERIMENTO BITANGENTE COPLANARE
            r_pi = a1*(1-e1); r_ai = a1*(1+e1);
            r_pf = a2*(1-e2); r_af = a2*(1+e2);
            
            % Gestione della geometria del trasferimento in base al type
            switch lower(maneuver)
                case 'pa'
                    r1 = r_pi; r2 = r_af;
                    th_start = 0; th_end = pi;
                case 'ap'
                    r1 = r_ai; r2 = r_pf;
                    th_start = pi; th_end = 2*pi;
                case 'pp'
                    % Controllo intersezione proiezioni XY
                    in1 = inpolygon(traj1(1,:), traj1(2,:), traj2(1,:), traj2(2,:));
                    in2 = inpolygon(traj2(1,:), traj2(2,:), traj1(1,:), traj1(2,:));
                    if ~(all(in1) || all(in2))
                        error('Manovra ''pp'' non consentita: le orbite proiettate sul piano xy si intersecano.');
                    end
                    
                    r1 = r_pi; r2 = r_pf;
                    th_start = 0; th_end = pi; 
                case 'aa'
                    % Controllo intersezione proiezioni XY
                    in1 = inpolygon(traj1(1,:), traj1(2,:), traj2(1,:), traj2(2,:));
                    in2 = inpolygon(traj2(1,:), traj2(2,:), traj1(1,:), traj1(2,:));
                    if (all(in1) || all(in2))
                        error('Manovra ''aa'' non consentita: le orbite proiettate sul piano xy NON si intersecano.');
                    end
                    
                    r1 = r_ai; r2 = r_af;
                    th_start = pi; th_end = 2*pi; 
            end
            
            a_t = (r1 + r2) / 2;
            e_t = abs(r2 - r1) / (r1 + r2);
            
            % Plot dell'arco di ellisse di trasferimento
            traj_t = get_traj(a_t, e_t, i1, OM1, om1, th_start, th_end);
            plot3(traj_t(1,:), traj_t(2,:), traj_t(3,:), 'r--', 'LineWidth', 2, 'DisplayName', ['Trasferimento Bitangente (' upper(maneuver) ')']);
            
        case 'bielliptic'
            % TRASFERIMENTO BIELLITTICO
            if isempty(varargin)
                error('Errore: Per il trasferimento biellittico devi fornire raggio apocentro (ra_t) come settimo argomento.');
            end
            ra_t = varargin{1};
            
            r_p1 = a1*(1-e1);
            r_p2 = a2*(1-e2);
            
            a_t1 = (r_p1 + ra_t) / 2;
            e_t1 = abs(ra_t - r_p1) / (r_p1 + ra_t);
            
            a_t2 = (r_p2 + ra_t) / 2;
            e_t2 = abs(ra_t - r_p2) / (r_p2 + ra_t);
            
            % Prima semi-ellisse (da pericentro iniziale ad apocentro ra_t)
            traj_t1 = get_traj(a_t1, e_t1, i1, OM1, om1, 0, pi);
            plot3(traj_t1(1,:), traj_t1(2,:), traj_t1(3,:), 'r--', 'LineWidth', 2, 'DisplayName', 'Trasferimento 1 (Biellittico)');
            
            % Seconda semi-ellisse (da apocentro ra_t a pericentro finale)
            traj_t2 = get_traj(a_t2, e_t2, i1, OM1, om1, pi, 2*pi);
            plot3(traj_t2(1,:), traj_t2(2,:), traj_t2(3,:), 'm--', 'LineWidth', 2, 'DisplayName', 'Trasferimento 2 (Biellittico)');
            
        case 'plane_change'
            % CAMBIO DI PIANO
            disp('Visualizzazione cambio di piano: le due orbite si intersecano ai nodi. Manovra istantanea.');
            
        case 'pericenter_change'
            % CAMBIO ANOMALIA PERICENTRO
            disp('Visualizzazione rotazione pericentro: le orbite mantengono stesso piano e forma. Manovra istantanea.');
            
        otherwise
            error('Manovra non riconosciuta. Usa: pa, ap, pp, aa, bielliptic, plane_change, pericenter_change.');
    end
    
    % --- 5. Disegna la Terra 3D
    R_earth = 6371; 
    [xE, yE, zE] = sphere(50); 
    
    try
        load topo topo topomap1; 
        surf(xE * R_earth, yE * R_earth, zE * R_earth, 'FaceColor', 'texturemap', ...
             'CData', topo, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        colormap(topomap1); 
    catch
        surf(xE * R_earth, yE * R_earth, zE * R_earth, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        colormap([0 0.2 0.6; 0.2 0.6 0.2; 0.6 0.4 0.2]); 
    end
    
    legend('show', 'Location', 'best');
    view(3);
    hold off;
    
    % =========================================================================
    % FUNZIONE ANNIDATA: Calcola i punti sfruttando la tua "par2car"
    % =========================================================================
    function rr_out = calc_traj(a, e, i, OM, om, th0, thf, dth, mu)
        th_vec = th0:dth:thf;
        if th_vec(end) ~= thf
            th_vec = [th_vec, thf];
        end
        rr_out = zeros(3, length(th_vec));
        for k = 1:length(th_vec)
            [rr_k, ~] = par2car(a, e, i, OM, om, th_vec(k), mu);
            rr_out(:, k) = rr_k;
        end
    end
end