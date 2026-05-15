function r_ijk = getPos(a, e, i, OM, om, th, mu)
   
    % function per animazione 3D in movimento nel codice
    % circolarizzaz_1impulso
        
    r_mag = (a*(1 - e^2)) / (1 + e*cos(th));
    r_pqw = [r_mag*cos(th); r_mag*sin(th); 0];
    
    R3_OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
    R1_i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
    r_ijk = R3_OM * R1_i * R3_om * r_pqw;
end