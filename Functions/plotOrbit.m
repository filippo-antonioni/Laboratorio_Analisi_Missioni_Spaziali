function plotOrbit(a, e, i, OM, om, th0, thf, dth, mu)
% 3D orbit plot
%
% plotOrbit(a, e, i, OM, om, th0, thf, dth, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% a           [1x1]   semi-major axis                  [km]
% e           [1x1]   eccentricity                     [-]
% i           [1x1]   inclination                      [rad]
% OM          [1x1]   RAAN                             [rad]
% om          [1x1]   pericenter anomaly               [rad]
% th0         [1x1]   initial true anomaly             [rad]
% thf         [1x1]   final true anomaly               [rad]
% dth         [1x1]   true anomaly step size           [rad]
% mu          [1x1]   gravitational parameter          [km^3/s^2]
th_vec=[th0:dth:thf];
% Controlla se l'ultimo elemento ha raggiunto esattamente thf. 
% Se non lo ha fatto, lo aggiungiamo a mano alla fine del vettore.
if th_vec(end) ~= thf
    th_vec = [th_vec, thf];
end
rr=[];
vv=[];
for k=1:length(th_vec)
[rr_i, vv_i] = par2car(a, e, i, OM, om, th_vec(k), mu);
rr=[rr rr_i];
vv=[vv vv_i];
end
    plot3(rr(1,:),rr(2,:),rr(3,:),'LineWidth',2)
    hold on;
    
    % --- AGGIUNTA: Disegna la Terra 3D con colori satellitari ---
    R_earth = 6371; % Raggio medio della Terra in km
    [xE, yE, zE] = sphere(50); % Crea una sfera con buona risoluzione (50x50)
    
    try
        % Carica la mappa topografica e i colori integrati in MATLAB
        load topo topo topomap1; 
        % Crea la superficie usando i dati topografici come texture
        surf(xE * R_earth, yE * R_earth, zE * R_earth, 'FaceColor', 'texturemap', ...
             'CData', topo, 'EdgeColor', 'none');
        colormap(topomap1); % Applica i colori: mari blu, terre verdi/marroni
    catch
        % Fallback base nel raro caso in cui i dati 'topo' manchino nell'installazione
        surf(xE * R_earth, yE * R_earth, zE * R_earth, 'EdgeColor', 'none');
        colormap([0 0.2 0.6; 0.2 0.6 0.2; 0.6 0.4 0.2]); % Palette base blu, verde, marrone
    end
    % --------------------------------------------------------------

    grid on;
    axis equal; 
    xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
end