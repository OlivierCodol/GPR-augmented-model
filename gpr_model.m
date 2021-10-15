function [HARD_FIT, SOFT_FIT, varargout] = gpr_model(W_D1_array, W_D2_array, Rw_array)



%-------------------------------
% GPR MODEL PARAMETERS
%-------------------------------
n_channel = 6;

W_SEL = 1;
W_CON = 1;
W_STN = 1;
W_SEL_GPi = -1;
W_CON_GPe = -1;
W_STN_GPi = 0.9;
W_STN_GPe = 0.9;
W_GPe_STN = -1;
W_GPe_GPi = -0.3;

e_SEL = 0.2;
e_CON = 0.2;
e_STN = -0.25;
e_GPe = -0.2;
e_GPi = -0.2;



%-------------------------------
% SIMULATION PARAMETERS
%-------------------------------
t = 2;                      % length of simulation
dt = 0.01;                  % time-step
n_steps = t / dt;

% selection threshold
th_SEL = 0;
th_DIST = 0.07;

% salience input time
ch1_onset = 0.6 / dt;
ch2_onset = 1.3 / dt;

% salience input value
salience1 = 0:.1:1; 
salience2 = 0:.1:1;
n_salience_1 = length(salience1);
n_salience_2 = length(salience2);

% selection times
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
if nargout > 2
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

for W_D1 = W_D1_array
    for W_D2 = W_D2_array
        
        % can use a parfor-loop below if desired to speed up simulations.
        % the code is already parfor-compatible.
        parfor k = 1:length(Rw_array) % parfor k = 1:length(Rw_array)
            
            % dopamine level
            Rw = Rw_array(k);
            LAMBDA = Rw2lambda(Rw);         

            % activity arrays
            a_SEL = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_CON = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_STN = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_GPe = zeros( n_channel , n_salience_1 , n_salience_2 );
            a_GPi = zeros( n_channel , n_salience_1 , n_salience_2 );
            
            % output arrays
            o_SEL = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_CON = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_STN = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_GPe = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            o_GPi = zeros( n_channel , n_salience_1 , n_salience_2 , n_steps );
            
            c = zeros( n_channel , n_salience_1 , n_salience_2 );
            
            % D1 and D2 scalers. This is to speed up the for-loop.
            D1_scaler = W_SEL * (1 + LAMBDA * W_D1);
            D2_scaler = W_CON * (1 - LAMBDA * W_D2);
            
            %%% SIMULATE MODEL
            for step = 2:n_steps
                % calculate salience changes
                if step == ch1_onset; c(1,:,:) = salience1_vec; end
                if step == ch2_onset; c(2,:,:) = salience2_vec; end
                
                % STRIATUM D1
                u_SEL = c * D1_scaler;
                a_SEL = (a_SEL - u_SEL) * decay_constant + u_SEL;
                o_SEL(:,:,:,step) = ramp(a_SEL,e_SEL,m);
                
                % STRIATUM D2
                u_CONT = c * D2_scaler;
                a_CON = (a_CON - u_CONT) * decay_constant + u_CONT;
                o_CON(:,:,:,step) = ramp(a_CON,e_CON,m);
                
                % STN
                u_STN = c * W_STN + W_GPe_STN * o_GPe(:,:,:,step-1);
                a_STN = (a_STN - u_STN) * decay_constant + u_STN;
                o_STN(:,:,:,step) = ramp(a_STN,e_STN,m);
                
                sum_o_STN = sum(o_STN(:,:,:,step));
                
                % GPe
                u_GPe =  sum_o_STN * W_STN_GPe + W_CON_GPe * o_CON(:,:,:,step);
                a_GPe = (a_GPe - u_GPe) * decay_constant + u_GPe;
                o_GPe(:,:,:,step) = ramp(a_GPe,e_GPe,m);
                
                % GPi
                u_GPi = sum_o_STN * W_STN_GPi + W_GPe_GPi * o_GPe(:,:,:,step) + W_SEL_GPi * o_SEL(:,:,:,step);
                a_GPi = (a_GPi - u_GPi) * decay_constant + u_GPi;
                o_GPi(:,:,:,step) = ramp(a_GPi,e_GPi,m);
            end
            
            % DETECT SIMULATION OUTCOME
            ch1_t_mid = o_GPi(1, :, :, t_mid);        % channel 1 t=1
            ch1_t_end = o_GPi(1, :, :, t_end);        % channel 1 t=2
            ch2_t_mid = o_GPi(2, :, :, t_mid);        % channel 2 t=1
            ch2_t_end = o_GPi(2, :, :, t_end);        % channel 2 t=2
            
            OUTCOME = 2*ones( n_salience_1 , n_salience_2 );
            OUTCOME( ch1_t_mid >th_SEL & ch1_t_end >th_SEL & ch2_t_end >th_SEL ) = 1;                       % no selection
            OUTCOME( ch1_t_mid<=th_SEL & ch1_t_end<=th_SEL & ch2_t_mid >th_SEL & ch2_t_end>th_DIST ) = 2;   % single channel selection (channel 1)
            OUTCOME( ch1_t_mid >th_SEL & ch1_t_end >th_SEL & ch2_t_end<=th_SEL & ch1_t_end>th_DIST ) = 2;   % single channel selection (channel 2)
            OUTCOME( ch1_t_mid<=th_SEL & ch1_t_end >th_SEL & ch2_t_end<=th_SEL & ch1_t_end>th_DIST ) = 2;   % channel switching
            OUTCOME( ch1_t_mid<=th_SEL & ch1_t_end >th_SEL & ch2_t_end >th_SEL ) = 1;                       % interference
            OUTCOME( ch1_t_mid<=th_SEL & ch1_t_end<=th_SEL & ch2_t_end<=th_SEL ) = 3;                       % dual channel selection
            OUTCOME( ch1_t_mid<=th_SEL & ch1_t_end<=th_SEL & ch2_t_end >th_SEL & ch2_t_end<=th_DIST ) = 2;  % distortion (channel 1 selected)
            OUTCOME( ch1_t_mid >th_SEL & ch1_t_end >th_SEL & ch2_t_end<=th_SEL & ch1_t_end<=th_DIST ) = 2;  % distortion (channel 2 selected)
            OUTCOME( ch1_t_mid<=th_SEL & ch1_t_end >th_SEL & ch2_t_end<=th_SEL & ch2_t_end<=th_DIST ) = 2;  % distortion (switching)
            
            if varargout_toggled
                OUTCOME_SAVE{k, iter} = OUTCOME;
            end
                
            % template fitting percentages
            SOFT_FIT(k, iter) = sum(sum(OUTCOME==OPTI_SOFT));
            HARD_FIT(k, iter) = sum(sum(OUTCOME==OPTI_HARD));
            
        end

        Progress = 100*iter/n_iter;
        fprintf('Progress ='); fprintf('%10.2f', Progress); fprintf(' '); fprintf('%%'); fprintf('\n');
        iter = iter + 1;
    end
end

% normalise & make into percentage
HARD_FIT = 100 * HARD_FIT / numel(OPTI_SOFT);
SOFT_FIT = 100 * SOFT_FIT / numel(OPTI_SOFT);

varargout = {OUTCOME_SAVE};





