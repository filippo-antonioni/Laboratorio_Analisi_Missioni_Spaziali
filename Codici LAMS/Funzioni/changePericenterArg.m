function [DeltaV, thi, thf] = changePericenterArg(a, e, omi, omf, mu)

% Change of Pericenter Argument maneuver
%
% [DeltaV, thetai, thetaf] = changePericenterArg(a, e, omi, omf, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% a             [1x1]   semi-major axis                     [km]
% e             [1x1]   eccentricity                        [-]
% omi           [1x1]   initial pericenter anomaly          [rad]
% omf           [1x1]   final pericenter anomaly            [rad]
% mu            [1x1]   gravitational parameter             [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% DeltaV        [1x1]   maneuver impulse                    [km/s]
% thi           [2x1]   initial true anomalies              [rad]
% thf           [2x1]   final true anomalies                [rad]
%

delta_om=omf-omi; %calcolo delta anomalia pericentro

thi_1=delta_om/2; % anomalia iniziale vera
thi_2=pi+delta_om/2;

thf_1=2*pi-delta_om/2; %anomalia finale vera
thf_2=pi-delta_om/2;

p=a*(1-e^2); %calcolo semilato retto
DeltaV=abs(2*sqrt(mu/p)*e*sin(delta_om/2)); %calcolo il deltaV per manovra
%non so se mettere abs per il deltav 
thi=mod([thi_1 thi_2]',2*pi);
thf=mod([thf_1 thf_2]',2*pi);

end