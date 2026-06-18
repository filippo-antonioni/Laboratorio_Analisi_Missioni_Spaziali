function [rr, vv] = par2car(a, e, i, OM, om, th, mu)

% Trasformation from Keplerian parameters to cartesian coordinates
%
% [rr, vv] = par2car(a, e, i, OM, om, th, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% a           [1x1]   semi-major axis                  [km]
% e           [1x1]   eccentricity                     [-]
% i           [1x1]   inclination                      [rad]
% OM          [1x1]   RAAN                             [rad]
% om          [1x1]   pericenter anomaly               [rad]
% th          [1x1]   true anomaly                     [rad]
% mu          [1x1]   gravitational parameter          [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% rr          [3x1]   position vector                  [km]
% vv          [3x1]   velocity vector                  [km/s]

p=a*(1-e^2); %semilato retto
%disp('CONTROLLA CHE GLI ANGOLI SIANO IN RADIANTI');

r = p / (1 + e * cos(th)); % distance from the central body

rr_pf = [r * cos(th); r * sin(th);0]; % position vector

vv_pf = sqrt(mu / p) * [-sin(th); e + cos(th); 0]; % velocity vector

% 1) Rotazione di OM (RAAN) intorno a k (asse Z)
R_OM = [ cos(OM),  sin(OM), 0;
        -sin(OM),  cos(OM), 0;
               0,        0, 1];

% 2) Rotazione di i (inclinazione) intorno a i' (asse X)
R_i = [1,       0,       0;
       0,  cos(i),  sin(i);
       0, -sin(i),  cos(i)];

% 3) Rotazione di om (anomalia del pericentro) intorno a k'' (asse Z)
R_om = [ cos(om),  sin(om), 0;
        -sin(om),  cos(om), 0;
               0,        0, 1];

T=R_om*R_i*R_OM; %matrice di rotazione totale da sistema eci a perifocale

rr=T'*rr_pf; %vettore posizione finale in eci

vv=T'*vv_pf; %vettore velocità finale in eci


end