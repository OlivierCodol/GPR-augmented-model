function [ratio, varargout] = compute_perf(measure, control)

ratio = max( measure ./ control, 0);% model perf vs control model perf ratio
tilde = log10( ratio );             % log10 of ratio
tilde( ratio==0 ) = nan;            % log10(x) for x=0 evaluates at -Inf
varargout = {tilde};

end