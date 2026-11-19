function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(ai, ei, af, ef, ra_t, mu)

% Bitangent transfer for elliptic orbits
%
% [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = bitangentTransfer(ai, ei, af, ef, type, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% ai            [1x1]   initial semi-major axis             [km]
% ei            [1x1]   initial eccentricity                [-]
% af            [1x1]   final semi-major axis               [km]
% ef            [1x1]   final eccentricity                  [-]
% ra_t          [1x1]   transfer orbits apocenter distance  [km]
% mu            [1x1]   gravitational parameter             [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV1       [1x1]   1st maneuver impulse                [km/s]
% DeltaV2       [1x1]   2nd maneuver impulse                [km/s]
% DeltaV3       [1x1]   3rd maneuver impulse                [km/s]
% Deltat1       [1x1]   maneuver time 1                     [s]
% Deltat2       [1x1]   maneuver time 2                     [s]
%
r_p1=ai*(1-ei); %raggio pericentro orbita trasferimento 1
r_p2=af*(1-ef); %raggio pericentro orbita trasferimento 2

a1=(r_p1+ra_t)/2; %semiasse orbita trasferimento 1
a2=(r_p2+ra_t)/2; %semiasse orbita trasferimento 2


DeltaV1=sqrt(mu)*(sqrt(2/r_p1-1/a1)-sqrt(2/r_p1-1/ai)); %deltav da orbita inziale a priam orbita trasferimento
DeltaV2=sqrt(mu)*(sqrt(2/ra_t-1/a2)-sqrt(2/ra_t-1/a1)); %deltav da orbita trasferimento 1 a orbita trasferimento 2
DeltaV3=sqrt(mu)*(sqrt(2/r_p2-1/af)-sqrt(2/r_p2-1/a2)); %deltav da orbita trasferimento 2 a orbita finale
%non so se mettere gli abs per i deltav, io direi di no visto che possono
%essere anche decelerazioni e quindi velocità negative, tuttavia se devo
%calcolare poi un deltav finale per stimare dei costi devo comunque
%sommarli con abs visto che inevitabilmente devo applicare delle variazioni
%di velocità, averle negative mi dà solo un senso fisico per capire se sto
%andando in una direzione o l'altra

Deltat1= pi*sqrt(a1^3/mu); %tempo trasferimento da orbita 1 a 2
Deltat2= pi*sqrt(a2^3/mu); %tempo trasferimento da 2 a finale



end