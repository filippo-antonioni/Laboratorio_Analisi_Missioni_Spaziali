function plotOrbit_circolarizzazione(a, e, i, OM, om, th_start, th_end, dth, mu, color, linestyle, linewidth, name)
    % plotOrbit_circolarizzazione: Traccia un arco di orbita
    th = th_start:dth:th_end;
    if th(end) ~= th_end
        th = [th, th_end];
    end

    r = (a*(1 - e^2)) ./ (1 + e*cos(th));

    % Perifocali
    r_pqw = [r.*cos(th); r.*sin(th); zeros(1, length(th))];

    % Matrici di rotazione
    R3_OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
    R1_i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];

    T_pqw2ijk = R3_OM * R1_i * R3_om;

    r_ijk = T_pqw2ijk * r_pqw;

    plot3(r_ijk(1,:), r_ijk(2,:), r_ijk(3,:), 'Color', color, 'LineStyle', linestyle, ...
          'LineWidth', linewidth, 'DisplayName', name);
end