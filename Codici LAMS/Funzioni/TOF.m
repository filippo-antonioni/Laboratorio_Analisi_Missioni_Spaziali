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
% th2         [1x1]  final true anomaly                           [rad]
% mu          [1x1]  gravitational parameter                      [km^3/s^2]
% 
% --------------------------------------------------------------------------
% Output arguments:
% deltat      [1x1]  time of flight                               [s]
% 
% -------------------------------------------------------------------------

th_vec = [th1 th2]'; 

% 1. Calcolo Anomalia Eccentrica
sin_E = (sqrt(1-e^2) .* sin(th_vec)) ./ (1 + e .* cos(th_vec));
cos_E = (e + cos(th_vec)) ./ (1 + e .* cos(th_vec));
E_vec = atan2(sin_E, cos_E);

% 2. Calcolo delle Anomalie Medie usando l'Eq. di Keplero (M = E - e*sin(E))
M_vec = E_vec - e .* sin(E_vec);

% 3. Variazione di Anomalia Media
delta_M = M_vec(2) - M_vec(1);

% 4. Modulo 2*pi sulla differenza
delta_M = mod(delta_M, 2*pi);

% 5. Calcolo del tempo di volo finale
deltat = delta_M * sqrt(a^3/mu);

end
