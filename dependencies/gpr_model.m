function [HARD_FIT, SOFT_FIT, varargout] = gpr_model(W_D1_array, W_D2_array, Rw_array, varargin)



%-------------------------------
% GPR MODEL PARAMETERS
%-------------------------------

% Synaptic Weights
weight_bounds = parsevarargin(varargin, 'weight_bounds', [1 1]);
W = get_weights(min(weight_bounds), max(weight_bounds));
W_SEL = W.W_SEL;
W_CON = W.W_CON;
W_CON_GPe_out = W.W_CON_GPe_out;
W_CON_GPe_inn = W.W_CON_GPe_inn;
W_CON_GPe_ta = W.W_CON_GPe_ta;
W_SEL_GPi = W.W_SEL_GPi;
W_STN = W.W_STN;
W_STN_GPe_out = W.W_STN_GPe_out;
W_STN_GPe_inn = W.W_STN_GPe_inn;
W_STN_GPe_ta = W.W_STN_GPe_ta;
W_STN_GPi = W.W_STN_GPi;
W_GPe_ta_CON = W.W_GPe_ta_CON;
W_GPe_ta_SEL = W.W_GPe_ta_SEL;
W_GPe_ta_GPe_ta = W.W_GPe_ta_GPe_ta;
W_GPe_out_GPe_out = W.W_GPe_out_GPe_out;
W_GPe_inn_GPe_inn = W.W_GPe_inn_GPe_inn;
W_GPe_out_GPe_inn = W.W_GPe_out_GPe_inn;
W_GPe_out_GPe_ta = W.W_GPe_out_GPe_ta;
W_GPe_inn_GPe_ta = W.W_GPe_inn_GPe_ta;
W_GPe_out_CON = W.W_GPe_out_CON;
W_GPe_out_SEL = W.W_GPe_out_SEL;
W_GPe_inn_CON = W.W_GPe_inn_CON;
W_GPe_inn_SEL = W.W_GPe_inn_SEL;
W_GPe_out_STN = W.W_GPe_out_STN;
W_GPe_out_GPi = W.W_GPe_out_GPi;
W_GPe_inn_STN = W.W_GPe_inn_STN;
W_GPe_inn_GPi = W.W_GPe_inn_GPi;

n_channel = 6;

e_SEL = 0.2;
e_CON = 0.2;
e_STN = -0.25;
e_GPe = -0.2;
e_GPi = -0.2;




%-------------------------------
% SIMULATION PARAMETERS
%-------------------------------
t = 0.6;                      % length of simulation
dt = 0.01;                  % time-step
n_steps = t / dt;

% selection threshold
th_SEL = 0;
th_DIST = 0.1 * 0.1032;  % 10% of tonic input value

% salience change time
ch2_onset = 0.3 / dt;

% salience input value
salience1 = 0:.1:1; 
salience2 = 0:.1:1;
n_salience_1 = length(salience1);
n_salience_2 = length(salience2);

% selection times (1 timestep before change of input)
t_mid = ch2_onset - 1;
t_end = n_steps - 1;

% activity and ramp output parameter
m = 1;                        % slope
decay_constant = exp(-25*dt); % 25 is gain 'k'



%-------------------------------
% INITIALIZE LOOP
%-------------------------------
[OPTI_HARD, OPTI_SOFT] = get_optimal_templates();

n_iter = length(W_D1_array)*length(W_D2_array);
iter = 1;

SOFT_FIT = nan(length(Rw_array), n_iter);
HARD_FIT = SOFT_FIT;
if nargout > 3
    OUTCOME_SAVE = cell(length(Rw_array), n_iter);
    varargout_toggled = true;
else
    OUTCOME_SAVE = {};
    varargout_toggled = false;
end

salience1_vec = repmat(salience1,[1,1,n_salience_2]);
salience2_vec = repmat(permute(salience2,[1,3,2]),[1,n_salience_1,1]);



%-------------------------------
% PERFORM SIMULATIONS
%-------------------------------
Progress = 100*iter/n_iter;
fprintf('Progress ='); fprintf('%10.2f', Progress); fprintf(' '); fprintf('%%');
        
for W_D1 = W_D1_array
    for W_D2 = W_D2_array
        
        % can use a parfor-loop below if desired to speed up simulations.
        % the code is already parfor-compatible.
        for k = 1:length(Rw_array) % parfor k = 1:length(Rw_array)
            
            % dopamine level
            Rw = Rw_array(k);
            LAMBDA = Rw2lambda(Rw);

            % activity arrays
            a_SEL = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_CON = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_STN = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_GPi = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_GPe_ta = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_GPe_out = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_GPe_inn = zeros( n_channel , n_salience_1 , n_salience_2 );
            
            
            % output arrays
            o_SEL = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_CON = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_STN = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_GPi = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_GPe_ta = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_GPe_out = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_GPe_inn = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            
            c = zeros( n_channel , n_salience_1 , n_salience_2 );
            c(1,:,:) = salience1_vec;
            
            % D1 and D2 scalers. This is to speed up the for-loop.
            D1_scaler = W_SEL * (1 + LAMBDA * W_D1);
            D2_scaler = W_CON * (1 - LAMBDA * W_D2);
            
            %%% SIMULATE MODEL
            for step = 2:n_steps
                % calculate salience changes
                %if step == ch1_onset; c(1,:,:) = salience1_vec; end
                if step == ch2_onset; c(2,:,:) = salience2_vec; end
                
                % STRIATUM D1
                u_SEL = c * D1_scaler + ...
                    W_GPe_ta_SEL .* sum(o_GPe_ta(:,:,:,step-1)) +...
                    W_GPe_out_SEL .* o_GPe_out(:,:,:,step-1) +...
                    W_GPe_inn_SEL .* o_GPe_inn(:,:,:,step-1);
                a_SEL = (a_SEL - u_SEL) * decay_constant + u_SEL;
                o_SEL(:,:,:,step) = ramp( a_SEL, e_SEL, m );
                
                % STRIATUM D2
                u_CON = c * D2_scaler + ...
                    W_GPe_ta_CON .* sum(o_GPe_ta(:,:,:,step-1)) + ...
                    W_GPe_out_CON .* o_GPe_out(:,:,:,step-1) + ...
                    W_GPe_inn_CON .* o_GPe_inn(:,:,:,step-1);
                a_CON = (a_CON - u_CON) * decay_constant + u_CON;
                o_CON(:,:,:,step) = ramp( a_CON, e_CON, m );
                
                % STN
                u_STN = c * W_STN + ...
                    W_GPe_out_STN .* o_GPe_out(:,:,:,step-1) +...
                    W_GPe_inn_STN .* o_GPe_inn(:,:,:,step-1);
                a_STN = (a_STN - u_STN) * decay_constant + u_STN;
                o_STN(:,:,:,step) = ramp( a_STN , e_STN, m );
                
                
                sum_o_STN = sum(o_STN(:,:,:,step));
                
                
                % GPe - Outer
                u_GPe_out = W_STN_GPe_out .* sum_o_STN + ...
                    W_CON_GPe_out .* o_CON(:,:,:,step) + ...
                    W_GPe_out_GPe_out .* o_GPe_out(:,:,:,step-1);
                a_GPe_out = (a_GPe_out - u_GPe_out) * decay_constant + u_GPe_out;
                o_GPe_out(:,:,:,step) = ramp( a_GPe_out, e_GPe, m );
                
                % GPe - Inner
                u_GPe_inn = W_STN_GPe_inn .* sum_o_STN + ...
                    W_CON_GPe_inn .* o_CON(:,:,:, step) + ...
                    W_GPe_out_GPe_inn .* o_GPe_out(:,:,:, step) + ...
                    W_GPe_inn_GPe_inn .* o_GPe_inn(:,:,:, step-1);
                a_GPe_inn = (a_GPe_inn - u_GPe_inn) * decay_constant + u_GPe_inn;
                o_GPe_inn(:,:,:,step) = ramp( a_GPe_inn, e_GPe, m );
                
                % GPe - TA
                u_GPe_ta = W_STN_GPe_ta .* sum_o_STN + ...
                    W_CON_GPe_ta .* o_CON(:,:,:, step) + ...
                    W_GPe_out_GPe_ta .* o_GPe_out(:,:,:, step) + ...
                    W_GPe_inn_GPe_ta .* o_GPe_inn(:,:,:, step) + ...
                    W_GPe_ta_GPe_ta .* o_GPe_ta(:,:,:, step-1);
                a_GPe_ta = (a_GPe_ta - u_GPe_ta) * decay_constant + u_GPe_ta;
                o_GPe_ta(:,:,:,step) = ramp( a_GPe_ta, e_GPe, m );
                

                % GPi
                u_GPi = W_STN_GPi .* sum_o_STN +...
                    W_SEL_GPi .* o_SEL(:,:,:,step) +...
                    W_GPe_out_GPi .* o_GPe_out(:,:,:,step) +...
                    W_GPe_inn_GPi .* o_GPe_inn(:,:,:,step);
                a_GPi = (a_GPi - u_GPi) * decay_constant + u_GPi;
                o_GPi(:,:,:,step) = ramp( a_GPi, e_GPi, m );
            end
            
            % DETECT SIMULATION OUTCOME
            ch1_mid = squeeze(o_GPi(1, :, :, t_mid));        % channel 1 t=1
            ch1_end = squeeze(o_GPi(1, :, :, t_end));        % channel 1 t=2
            ch2_end = squeeze(o_GPi(2, :, :, t_end));        % channel 2 t=2
            
            OUTCOME = ones( n_salience_1 , n_salience_2 );

            OUTCOME( ch1_mid >  th_SEL & ch1_end  > th_SEL & ch2_end  > th_SEL ) = 1;                         % no selection
            OUTCOME( ch1_mid <= th_SEL & ch1_end <= th_SEL & ch2_end  > th_SEL & ch2_end > th_DIST ) = 2;   % single channel selection (channel 1)
            OUTCOME( ch1_mid >  th_SEL & ch1_end  > th_SEL & ch2_end <= th_SEL & ch1_end > th_DIST ) = 2;   % single channel selection (channel 2)
            OUTCOME( ch1_mid <= th_SEL & ch1_end  > th_SEL & ch2_end <= th_SEL & ch1_end > th_DIST ) = 3;   % channel switching
            OUTCOME( ch1_mid <= th_SEL & ch1_end  > th_SEL & ch2_end  > th_SEL ) = 4;                         % interference
            OUTCOME( ch1_mid <= th_SEL & ch1_end <= th_SEL & ch2_end <= th_SEL ) = 5;                         % dual channel selection
            OUTCOME( ch1_mid <= th_SEL & ch1_end <= th_SEL & ch2_end  > th_SEL & ch2_end <= th_DIST ) = 6;  % distortion (channel 1 selected)
            OUTCOME( ch1_mid >  th_SEL & ch1_end  > th_SEL & ch2_end <= th_SEL & ch1_end <= th_DIST ) = 6;  % distortion (channel 2 selected)
            OUTCOME( ch1_mid <= th_SEL & ch1_end  > th_SEL & ch2_end <= th_SEL & ch1_end <= th_DIST ) = 6;  % distortion (switching)

            if varargout_toggled
                OUTCOME_SAVE{k, iter} = OUTCOME;
            end

            % template fitting percentages
            SOFT_FIT(k, iter) = sum(sum(OUTCOME==OPTI_SOFT));
            HARD_FIT(k, iter) = sum(sum(OUTCOME==OPTI_HARD));
            
        end

        Progress = 100*iter/n_iter;
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b')
        fprintf('%10.2f', Progress); fprintf(' '); fprintf('%%');
        iter = iter + 1;
    end
end

% normalise & make into percentage
HARD_FIT = 100 * HARD_FIT / numel(OPTI_SOFT);
SOFT_FIT = 100 * SOFT_FIT / numel(OPTI_SOFT);

varargout = {W, OUTCOME_SAVE};
fprintf('\n');




