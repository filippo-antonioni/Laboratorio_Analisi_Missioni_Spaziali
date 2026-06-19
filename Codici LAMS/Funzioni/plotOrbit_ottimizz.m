function plotOrbit_ottimizz(a, e, i, OM, om, mu, th_vec, col, stile, width, nome)
    pts = zeros(3, length(th_vec));
    for j = 1:length(th_vec)
        [pts(:,j), ~] = par2car(a, e, i, OM, om, th_vec(j), mu);
    end
    plot3(pts(1,:), pts(2,:), pts(3,:), 'Color', col, 'LineStyle', stile, ...
          'LineWidth', width, 'DisplayName', ['\textbf{', nome, '}']);
end