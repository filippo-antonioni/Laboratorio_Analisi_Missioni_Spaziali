function stop = outfun_fmincon(x, optimValues, state)
    global history_fmincon
    stop = false;
    if strcmp(state, 'iter') || strcmp(state, 'init')
        history_fmincon = [history_fmincon; x(1), x(2), x(3), optimValues.fval];
    end
end