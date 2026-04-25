function stop = outfun_pso(optimValues, state)
    global history_pso
    stop = false;
    % Salviamo: [Iterazione, Miglior Valore, Valori di tutte le particelle]
    history_pso = [history_pso; optimValues.iteration, optimValues.bestfval, optimValues.swarmfvals'];
end