function plotOrbit_sole(a, e, i, OM, om, th0, thf, dth, mu)
% 3D orbit plot
%
% plotOrbit_sole(a, e, i, OM, om, th0, thf, dth, mu)
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

th_vec = [th0:dth:thf];

% Controlla se l'ultimo elemento ha raggiunto esattamente thf. 
% Se non lo ha fatto, lo aggiungiamo a mano alla fine del vettore.
if th_vec(end) ~= thf
    th_vec = [th_vec, thf];
end

rr = [];
vv = [];
for k = 1:length(th_vec)
    [rr_i, vv_i] = par2car(a, e, i, OM, om, th_vec(k), mu);
    rr = [rr, rr_i];
    vv = [vv, vv_i];
end

plot3(rr(1,:), rr(2,:), rr(3,:), 'LineWidth', 2)
hold on;

% --- AGGIUNTA: Disegna il Sole 3D ---
R_sun = 696340; % Raggio medio del Sole in km
[xS, yS, zS] = sphere(50); % Crea una sfera con buona risoluzione (50x50)

% Crea la superficie usando un colore giallo/arancione
surf(xS * R_sun, yS * R_sun, zS * R_sun, 'FaceColor', [1, 0.8, 0], 'EdgeColor', 'none');

% Aggiunge l'illuminazione per dare un volume 3D alla sfera
camlight;
lighting gouraud; 
% --------------------------------------------------------------

grid on;
axis equal; 
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

end