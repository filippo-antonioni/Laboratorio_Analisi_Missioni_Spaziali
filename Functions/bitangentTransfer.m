function [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)

% Bitangent transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(ai, ei, af, ef, type, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai            [1x1]   initial semi-major axis             [km]
% ei            [1x1]   initial eccentricity                [-]
% af            [1x1]   final semi-major axis               [km]
% ef            [1x1]   final eccentricity                  [-]
% type          [char]  maneuver type
% mu            [1x1]   gravitational parameter             [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV1       [1x1]   1st maneuver impulse                [km/s]
% DeltaV2       [1x1]   2nd maneuver impulse                [km/s]
% Deltat        [1x1]   maneuver time                       [s]
%

r_pi=a_i*(1-e_i); %raggio pericentro orbita iniziale
r_af=a_f*(1+e_f); %raggio apocentro orbita finale
r_ai=a_i*(1+e_i); %raggio apocentro orbita iniziale
r_pf=a_f*(1-e_f); %raggio pericentro orbita finale

switch type
    case 'pa'
        
        r1=r_pi;
        r2=r_af;
        
    case 'ap'
        r1=r_ai;
        r2=r_pf;
        
    case 'pp'
        r1=r_pi;
        r2=r_pf;
    case 'aa'
        r1=r_ai;
        r2=r_af;
    otherwise
        error('Hai sbagliato a inserire il type...');
end
    a=(r1+r2)/2; %calcolo semiasse maggiore
    DeltaV1=sqrt(mu)*(sqrt(2/r1-1/a)-sqrt(2/r1-1/a_i)); %calcolo il primo delta_V
    
    DeltaV2=sqrt(mu)*(sqrt(2/r2-1/a_f)-sqrt(2/r2-1/a)); %calcolo il secondo delta_V
    
    Deltat=pi*sqrt(a^3/mu); %calcolo il tempo di manovra

end