function [Wx, d_Fh, d_Fs, H_max, S_max ] = compute_proxies(hard_array , soft_array, Rw_array)

n_combinations = size(hard_array,2);
[Wx, d_Fh, d_Fs, H_max, S_max ] = deal( zeros(n_combinations, 1) );


% This assumes that soft & hard Template Fitting curves intersect only once
for combi = 1:n_combinations
    
    %------------------------
    % locate max template fits
    %------------------------
    hard = hard_array(:,combi);
    soft = soft_array(:,combi);
    H_max(combi) = max(hard);
    S_max(combi) = max(soft);
    
    %------------------------
    % compute integrals (trapz method)
    %------------------------
    d_fit = hard - soft;
    d_Rw = diff(Rw_array(:));
    I = 0.5 * d_Rw .* (d_fit(1:end-1) + d_fit(2:end));
    d_Fh(combi) =     sum( I(I>0) )  ./ sum( d_Rw(I>0) );
    d_Fs(combi) = abs(sum( I(I<0) )) ./ sum( d_Rw(I<0) );
    
    %------------------------
    % find crosspoint Wx
    %------------------------
    sign_change = diff(sign(d_fit));
    takeoff = find( sign_change~=0 , 1);  % first time a difference comes up
    
    if isempty(takeoff)  % hard and soft are always equal to each other, no crosspoint
        Wx(combi) = nan;
    elseif sign_change(takeoff) < 0  % once a difference arise, we are immediately in soft regime
        Wx(combi) = Rw_array(takeoff + 1);
    else % we start with hard regime
        crosspoint = find( sign_change(takeoff:end) < 0, 1, 'last'); % find the last time we switch to soft regime
        
        if isempty(crosspoint) % if we never switch, no crosspoint
            Wx(combi) = nan;
        else
            Wx(combi) = Rw_array(crosspoint + takeoff - 1);
        end
    end

end


