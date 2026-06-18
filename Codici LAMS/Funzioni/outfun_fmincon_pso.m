function stop = outfun_fmincon_pso(x, optimValues, state)
    global history_fmincon
    stop = false;
    if strcmp(state, 'iter') || strcmp(state, 'init')
        history_fmincon = [history_fmincon; optimValues.iteration, optimValues.fval];
    end
end