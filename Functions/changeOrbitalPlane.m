function [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)

% Change of Plane maneuver
%
% [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% a         [1x1]   semi-major axis                 [km]
% e         [1x1]   eccentricity                    [-]
% i_i       [1x1]   initial inclination             [rad]
% OMi       [1x1]   initial RAAN                    [rad]
% omi       [1x1]   initial pericenter anomaly      [rad]
% i_f       [1x1]   final inclination               [rad]
% OMf       [1x1]   final RAAN                      [rad]
% mu        [1x1]   gravitational parameter         [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV    [1x1]   maneuver impulse                [km/s]
% omf       [1x1]   final pericenter anomaly        [rad]
% theta     [1x1]   true anomaly at maneuver        [rad]
% 

    % 1. Calcolo delle variazioni
    dOM = OMf - OMi;
    di  = i_f - i_i;
    
    % 2. Determinazione del caso
    if dOM > 0 && di > 0
        caso = 1;
    elseif dOM > 0 && di < 0
        caso = 2;
    elseif dOM < 0 && di > 0
        caso = 3;
    elseif dOM < 0 && di < 0
        caso = 4;
    else
        error('Variazioni nulle (dOM=0 o di=0) non gestite in questo algoritmo.');
    end
    
    % 3. Calcoli comuni a tutti i casi (Triangolo sferico)
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(dOM));
    
    abs_dOM = abs(dOM);
    sin_ui = (sin(abs_dOM) / sin(alpha)) * sin(i_f);
    sin_uf = (sin(abs_dOM) / sin(alpha)) * sin(i_i);
    
    % 4. Switch per le logiche specifiche dei 4 casi
    switch caso
        case 1
            cos_ui =  (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
            cos_uf =  (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
            
            u_i = atan2(sin_ui, cos_ui);
            u_f = atan2(sin_uf, cos_uf);
            
            theta = u_i - omi;
            omf   = u_f - theta;
            
        case 2
            cos_ui =  (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
            cos_uf = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
            
            u_i = atan2(sin_ui, cos_ui);
            u_f = atan2(sin_uf, cos_uf);
            
            theta = 2*pi - u_i - omi;
            omf   = 2*pi - u_f - theta;
            
        case 3
            cos_ui = (-cos(i_f) + cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
            cos_uf =  (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
            
            u_i = atan2(sin_ui, cos_ui);
            u_f = atan2(sin_uf, cos_uf);
            
            theta = 2*pi - u_i - omi;
            omf   = 2*pi - u_f - theta;
            
        case 4
            cos_ui =  (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
            cos_uf = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
            
            u_i = atan2(sin_ui, cos_ui);
            u_f = atan2(sin_uf, cos_uf);
            
            theta = u_i - omi;
            omf   = u_f - theta;
    end
    
    % Mantenere gli angoli nel range [0, 2*pi]
    theta = mod(theta, 2*pi);
    omf   = mod(omf, 2*pi);
    if cos(theta)>=0
        theta= theta + pi;
    end

    % 5. Calcolo della DeltaV
    % Calcolo del semi-lato retto
    p = a * (1 - e^2);
    
    % Componente trasversale della velocità al nodo
    v_theta = sqrt(mu/p)*(1+e*cos(theta)); 
    
    % Nel cambio di piano puro (a ed e costanti), ruota solo v_theta
    DeltaV = 2 * v_theta * sin(alpha / 2);
           
end