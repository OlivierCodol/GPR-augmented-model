function [Wx, d_Fh, d_Fs, H_max, S_max ] = compute_proxies(hard_array , soft_array, Rw_array)

n_combinations = size(hard_array,2);
[Wx, d_Fh, d_Fs, H_max, S_max ] = deal( zeros(n_combinations, 1) );


% This assumes that soft & hard Template Fitting curves intersect only once
for combi = 1:n_combinations
    
    hard = hard_array(:,combi);
    soft = soft_array(:,combi);
    fit_diff = hard - soft;
    
    H_max(combi) = max(hard);
    S_max(combi) = max(soft);
    
    crosspoints = find( diff(sign(fit_diff))~=0 ); % find all crossing points

    if ~isempty(crosspoints)
        cross_first = crosspoints(1);
        cross_last = crosspoints(end);
        Wx(combi) = Rw_array ( cross_first );
        
        x = Rw_array(1:cross_first+1);
        y = [fit_diff(1:cross_first) ; 0]; % defensive padding against 1-element vectors
        d_Fh(combi) = trapz(x,y) / (x(end)-x(1));
        
        x = Rw_array(cross_last:end);
        y = -fit_diff(cross_last:end);
        d_Fs(combi) = trapz(x,y) / (x(end)-x(1));
    else
        x = Rw_array;
        y = fit_diff;
        if fit_diff(1)>0
            d_Fh(combi) = trapz(x,y) / (x(end)-x(1));
        else
            d_Fs(combi) = trapz(x,y) / (x(end)-x(1));
        end
    end
	
end

