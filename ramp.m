%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DATE: 24/10/2020
%%%% WHAT: piece-wise linear output function (i.e. ramp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function output = ramp(a,e,m)

%---- for 'a' a scalar input ---
% if a < e
%     output = 0;
% elseif a <= 1/m + e
%     output = m * (a - e);
% else 
%     output = 1;
% end

%---- for 'a' an array of values ---
if m==1 % optimised for speed
    output = a - e;
    output = min(output,1);
    output = max(output,0);
else
    case_zero = (a < e);
    case_a = (a >= e & a <= 1/m +e);
    case_one = (a > 1/m + e);
    
    output = nan(size(a));
    output(case_zero) = 0;
    output(case_a) = m .* (a(case_a) - e);
    output(case_one) = 1;
end
