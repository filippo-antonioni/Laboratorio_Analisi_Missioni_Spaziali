function [a, e, i, OM, om, th] = car2par(rr, vv, mu)

% Trasformation from cartesian coordinates to Keplerian parameters
%
% [a, e, i, OM, om, th] = car2par(rr, vv, mu)
%
% -------------------------------------------------------------------------
% Input arguments:
% rr          [3x1]   position vector                  [km]
% vv          [3x1]   velocity vector                  [km/s]
% mu          [1x1]   gravitational parameter          [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% a           [1x1]   semi-major axis                  [km]
% e           [1x1]   eccentricity                     [-]
% i           [1x1]   inclination                      [rad]
% OM          [1x1]   RAAN                             [rad]
% om          [1x1]   pericenter anomaly               [rad]
% th          [1x1]   true anomaly                     [rad]
%
% -------------------------------------------------------------------------

r=norm(rr); %calcolo il modulo del vettore posizione
v = norm(vv); %calcolo il modulo del vettore velocità

epsilon= 0.5*v^2-mu/r; %energia meccanica specfica  

a=(2/r-v^2/mu)^-1; %calcolo il semiasse maggiore

hh=cross(rr,vv); %calcolo vettore momento angolare specifico
h=norm(hh); %calcolo il modulo

ee=cross(vv,hh)/mu-rr/r; %calcolo vettore eccentricità
e=norm(ee); %calcolo eccentricità

kk=[0 0 1]';
i=acos(hh(3)/h); %calcolo inclinazione orbita, prendo terza componente del vettore hh

NN=cross(kk,hh)/norm(cross(kk,hh)); %Linea dei nodi

if NN(2)>=0

    OM=acos(NN(1)); %calcolo RAAN nel caso in cui angolo 0<OM<pi
    
else 
    OM=2*pi-acos(NN(1)); %calcolo RAAN nel caso in cui angolo sia pi<OM<2pi

end

if ee(3)>=0
    om=acos(dot(NN,ee)/e); %calcolo anomalia del pericentro caso 0<om<pi

else 
    om=2*pi-acos(dot(NN,ee)/e); %calcolo anomalia del pericentro caso pi<om<2pi

end

v_r=dot(vv,rr)/r; %velocità radiale in modulo
 
 if v_r>=0 
    
     th=acos(dot(rr,ee)/(r*e)); %anomalia vera caso in  cui ci allontaniamo da pericentro 0<theta<pi

 else 
    th=2*pi-acos(dot(rr,ee)/(r*e)); %anomalia vera avvicinandosi pi<theta<2pi
 
 end



end