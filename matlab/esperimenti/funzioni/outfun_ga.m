function [state, options, optchanged] = outfun_ga(options, state, flag)
    global history_ga
    optchanged = false;
    if strcmp(flag, 'iter') || strcmp(flag, 'init')
        history_ga = [history_ga; state.Population, state.Score];
    end
end