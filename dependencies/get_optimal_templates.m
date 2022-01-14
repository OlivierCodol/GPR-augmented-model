function [OPTI_HARD, OPTI_SOFT] = get_optimal_templates()

OPTI_HARD = ones(11, 11);
OPTI_HARD(1:8, :) = 2;
OPTI_HARD(9:11, 4:11) = 2;
for j = 1:8
    OPTI_HARD (j, 12-j) = 4;
    OPTI_HARD (j, 3:11-j) = 3;
end
% swap salience 1 and 2 and flip axis to match matrix dimension order
OPTI_HARD = permute(flipud(OPTI_HARD), [2,1]);

OPTI_SOFT = ones(11, 11);
OPTI_SOFT(1:8, :) = 2;
OPTI_SOFT(:, 4:11) = 2;
OPTI_SOFT(1:8, 4:11) = 5;
OPTI_SOFT = permute(flipud(OPTI_SOFT), [2,1]);
