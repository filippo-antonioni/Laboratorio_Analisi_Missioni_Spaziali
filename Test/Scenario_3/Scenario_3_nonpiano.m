clear
clc
close all

%% SECTION 1 - Geocentrica -> Eliocentrica (Iperbole di Fuga)
%%% Load dati orbita finale geocentrica (scenario 1) - ORBITA DI PARCHEGGIO %%%
T = load("DatiSC1-2026.txt");
mu_earth = 398600;
m_earth = 5.974e24; % [kg]
gruppo = 23;
R_earth=6371;

% Orbita target/finale (vettori di stato)
rr_f = T(gruppo,8:10)';
vv_f = T(gruppo,11:13)';
state_f = [rr_f; vv_f];

% Trasformazione per avere i parametri kepleriani bersaglio (corretto: mu_earth)
[a_op, e_op, i_op, OM_op, om_op, th_op] = car2par(rr_f, vv_f, mu_earth);

% Moduli del raggio e velocità di pericentro dell'orbita di parcheggio
r_op_p = a_op * (1 - e_op);
v_op = sqrt(mu_earth * (2/r_op_p - 1/a_op));

% Modulo del raggio di pericentro dell'orbita di parcheggio
r_op_a = a_op * (1 + e_op);

%%% Load dati orbita eliocentrica (scenario 2) %%%
load('workspace_sc2_gridsearch.mat'); 
mu_sun  = 1.32712440018e11; % Parametro gravitazionale del Sole [km^3/s^2]
m_sun   = 1.989e30;         % [kg]
R_sun   = 696340;           % Raggio del Sole [km]
AU      = 149597870.7;      % 1 Unità Astronomica [km]

%%% Calcolo SOI terra (usando la distanza reale alla partenza) %%%
% r1_opt è stato calcolato nello Scenario 2 (o lo ricalcoli con ottimo_th1)
R_terra_partenza = norm(r1_opt); 
r_SOI_terra = R_terra_partenza * (m_earth/m_sun)^(2/5);

% Velocità eliocentrica della Terra alla partenza (V1i)
[r1i, V1i] = par2car(a_T, e_T, i_T, OM_T, om_T, ottimo_th1, mu_sun);

% Velocità eliocentrica del satellite sull'orbita di trasferimento alla partenza (V1T)
[R1T, V1T] = par2car(ottimo_aT, ottimo_eT, i_trasf_opt, OM_transf_opt, ottimo_omT, th1_T_opt, mu_sun);

% Eccesso iperbolico vettoriale e modulo
v_inf_1_vec = V1T - V1i;
%v_inf_1 = norm(v_inf_1_vec);


% Matrice di rotazione da SDR eclittico a SDR ECI
epsilon = deg2rad(23.45);
T = [1 0 0;
   0 cos(epsilon) sin(epsilon);
   0 -sin(epsilon) cos(epsilon)];

% Converto la v_inf da eclittico a ECI 
v_inf_1_vec=T'*v_inf_1_vec;
v_inf_1 = norm(v_inf_1_vec);

% Trovo versore che definisce l'asintoto uscente dell'iperbole 

r_inf_vers = v_inf_1_vec / v_inf_1;

% Procedo a calcolare con un for il raggio r_h, assumo come punto di
% partenza il raggio di pericentro dell'orbita di parcheggio essendo il
% punto sull'orbita di parcheggio a distanza minima dal fuoco

th_h_vec=linspace(0,2*pi,1000);
N_r = length(th_h_vec);

% Calcolo semiasse maggiore iperbole
a_h=-mu_earth/v_inf_1^2;

% Vettori che uso nel ciclo for
dv_1_vec=zeros(N_r,1);
r_h_vec=zeros(N_r,3);
e_h_vec=zeros(N_r,1);
th_h_vec=zeros(N_r,1);
th_inf_vec=zeros(N_r,1);


for i=1:N_r 
    
    % Calcolo r_h_i associato al th_h_i
    [r_h_i, v_h_i] = par2car(a_op, e_op, i_op, OM_op, om_op, th_h_vec(i), mu_earth);
    
    % Definisco vettore per risolvere il sistema con Newton
    % x=[e_h,th_h,th_inf]

    r_h_i_norm=norm(r_h_i);
    r_h_i_vers= r_h_i / r_h_i_norm;
    
    % Calcolo e_h guess iniziale sufficientemente vicino spero
    
    e_h_0= 1-r_h_i_norm/a_h;
    th_inf_0 = acos(-1 / e_h_0);
    th_h_0 = th_inf_0 - alfa;

    
    x0=[e_h_0;th_inf_0;th_h_0];
    
    alfa=acos(dot(r_h_i_vers,r_inf_vers));
    
    fun=@(x) [r_h_i_norm - a_h*(1-x(1)^2)/(1+x(1)*cos(x(2)));
             cos(x(3)) + 1/x(1);
             x(3) - x(2) -alfa];   
    
    x_i=fsolve(fun,x0);
    
    % Calcolo raggio di pericentro iperbole i-esima che ho trovato
    r_p_h_i= a_h*(1-x_i(1));

    % Controllo vincoli
    if (x_i(3)>pi/2 && x_i(3)<pi) && r_p_h_i>R_earth
        
        r_h_vec(i)=r_h_i;
        e_h_vec(i)=x_i(1);
        th_h_vec(i)=x_i(2);
        th_inf_vec=x_i(3);
        

    end 




end


