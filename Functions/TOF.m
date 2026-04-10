function deltat=TOF(a,e,th1,th2,mu)
% Time of flight 
% 
% deltat=TOF(a,e,th1,th2,mu)
% 
% --------------------------------------------------------------------------
% Input arguments:
% a           [1x1]  semi-major axis                              [km]
% e           [1x1]  eccentricity                                 [-]
% th1         [1x1]  initial true anomaly                         [rad]  
% th2         [1x1]  final trur anomaly                           [rad]
% mu          [1x1]  gravitational parameter                      [km^3/s^2]
% 
% --------------------------------------------------------------------------
% Output arguments:
% deltat      [1x1]  time of flight                               [s]
% 
% -------------------------------------------------------------------------
th_vec = [th1 th2]'; % vettore anomalie vere 

% --- MODIFICA: Uso di atan2 per evitare salti di quadrante ---
sin_E = (sqrt(1-e^2) .* sin(th_vec)) ./ (1 + e .* cos(th_vec));
cos_E = (e + cos(th_vec)) ./ (1 + e .* cos(th_vec));
E_vec = atan2(sin_E, cos_E);

% Riportiamo le anomalie eccentriche nel range positivo [0, 2*pi]
E_vec = mod(E_vec, 2*pi);
% -------------------------------------------------------------

% Calcolo del tempo di volo tramite Equazione di Keplero
deltat = sqrt(a^3/mu)*(E_vec(2) - E_vec(1) - e*(sin(E_vec(2)) - sin(E_vec(1))));

if th1 > th2
    T = 2*pi*sqrt(a^3/mu);
    deltat = deltat + T;
end
end