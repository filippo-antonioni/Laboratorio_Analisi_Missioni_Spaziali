function plotOrbit_circolarizzazioneStaticAnim(a, e, i, OM, om, mu, col, stile, width, nome)
    
        % function per circolarizzazione per plot dell'animazione 


    th = linspace(0, 2*pi, 200);
    pts = zeros(3, 200);
    for j = 1:200
        r_mag = (a*(1 - e^2)) / (1 + e*cos(th(j)));
        rpqw = [r_mag*cos(th(j)); r_mag*sin(th(j)); 0];
        R3OM = [cos(OM) -sin(OM) 0; sin(OM) cos(OM) 0; 0 0 1];
        R1i  = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
        R3om = [cos(om) -sin(om) 0; sin(om) cos(om) 0; 0 0 1];
        pts(:,j) = R3OM * R1i * R3om * rpqw;
    end
    plot3(pts(1,:), pts(2,:), pts(3,:), 'Color', [col 0.4], 'LineStyle', stile, 'LineWidth', width, 'DisplayName', nome);
end