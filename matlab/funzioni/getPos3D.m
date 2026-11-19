function r_ijk = getPos3D(a, e, i, OM, om, th)
    r = (a*(1 - e^2)) / (1 + e*cos(th));
    
    % Vettore in componenti perifocali (piano dell'orbita)
    r_pqw = [r*cos(th); r*sin(th); 0];
    
    % Matrici di rotazione
    R3_OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
    R1_i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
    
    % Rotazione da perifocale a inerziale geocentrico (ECI)
    T_pqw2ijk = R3_OM * R1_i * R3_om;
    r_ijk = T_pqw2ijk * r_pqw;
end
