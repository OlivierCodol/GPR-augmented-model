function [c, o_GPi] = dynamics_sim(W_D1, W_D2, Rw)

%-------------------------------
% GPR MODEL PARAMETERS
%-------------------------------
n_channel = 6;





e_SEL = 0.2;
e_CON = 0.2;
e_STN = -0.25;
e_GPe = -0.2;
e_GPi = -0.2;

% Synaptic Weights
W_SEL = 1;
W_CON = 1;
W_CON_GPe_out = -.9;
W_CON_GPe_inn = -.9;
W_CON_GPe_ta = -.9;
W_SEL_GPi = -1;
W_STN = 1;
W_STN_GPe_out = 0.8;
W_STN_GPe_inn = 0.8;
W_STN_GPe_ta = 0.8;
W_STN_GPi = 0.9;
W_GPe_ta_CON = -0.25;
W_GPe_ta_SEL = -0.25;
W_GPe_ta_GPe_ta = -0.75;
W_GPe_out_GPe_out = -0.75;
W_GPe_inn_GPe_inn = -0.75;
W_GPe_out_GPe_inn = -0.3;
W_GPe_out_GPe_ta = -0.75;
W_GPe_inn_GPe_ta = -0.75;
W_GPe_out_CON = 0.5;
W_GPe_out_SEL = 0.5;
W_GPe_inn_CON = 0.25;
W_GPe_inn_SEL = 0.25;
W_GPe_out_STN = -0.8;
W_GPe_out_GPi = -1;
W_GPe_inn_STN = -0.8;
W_GPe_inn_GPi = -0.2;



%-------------------------------
% SIMULATION PARAMETERS
%-------------------------------
t = 3;                      % length of simulation
dt = 0.01;                  % time-step
n_steps = t / dt;

% salience input time
ch1_onset = 1 / dt;
ch2_onset = 2 / dt;

% activity and ramp output parameter
m = 1;                        % slope
decay_constant = exp(-25*dt); % 25 is gain 'k'



%-------------------------------
% PERFORM SIMULATIONS
%-------------------------------

% dopamine level
LAMBDA = Rw2lambda(Rw);

% activity arrays
a_SEL = zeros( n_channel , 1 );
a_CON = zeros( n_channel , 1 );
a_STN = zeros( n_channel , 1 );
a_GPi = zeros( n_channel , 1 );
a_GPe_ta = zeros( n_channel , 1 );
a_GPe_out = zeros( n_channel , 1 );
a_GPe_inn = zeros( n_channel , 1 );

% output arrays
o_SEL = zeros( n_channel , n_steps );
o_CON = zeros( n_channel , n_steps );
o_STN = zeros( n_channel , n_steps );
o_GPi = zeros( n_channel , n_steps );
o_GPe_ta = zeros( n_channel , n_steps );
o_GPe_out = zeros( n_channel , n_steps );
o_GPe_inn = zeros( n_channel , n_steps );

c = zeros( n_channel , n_steps );

% D1 and D2 scalers. This is to speed up the for-loop.
D1_scaler = W_SEL * (1 + LAMBDA * W_D1);
D2_scaler = W_CON * (1 - LAMBDA * W_D2);

%%% SIMULATE MODEL
for step = 2:n_steps
    % calculate salience changes
    if step == ch1_onset; c(1, step:end) = 0.4; end
    if step == ch2_onset; c(2, step:end) = 0.6; end
    
    % STRIATUM D1
    u_SEL = c(:, step) * D1_scaler + ...
        W_GPe_ta_SEL .* sum(o_GPe_ta(:,step-1)) +...
        W_GPe_out_SEL .* o_GPe_out(:,step-1) +...
        W_GPe_inn_SEL .* o_GPe_inn(:,step-1);
    a_SEL = (a_SEL - u_SEL) * decay_constant + u_SEL;
    o_SEL(:,step) = ramp( a_SEL, e_SEL, m );
    
    % STRIATUM D2
    u_CON = c(:, step) * D2_scaler + ...
        W_GPe_ta_CON .* sum(o_GPe_ta(:,step-1)) + ...
        W_GPe_out_CON .* o_GPe_out(:,step-1) + ...
        W_GPe_inn_CON .* o_GPe_inn(:,step-1);
    a_CON = (a_CON - u_CON) * decay_constant + u_CON;
    o_CON(:,step) = ramp( a_CON, e_CON, m );
    
    % STN
    u_STN = c(:, step) * W_STN + ...
        W_GPe_out_STN .* o_GPe_out(:,step-1) +...
        W_GPe_inn_STN .* o_GPe_inn(:,step-1);
    a_STN = (a_STN - u_STN) * decay_constant + u_STN;
    o_STN(:,step) = ramp( a_STN , e_STN, m );
    
    
    sum_o_STN = sum(o_STN(:,step));
    
    
    % GPe - Outer
    u_GPe_out = W_STN_GPe_out .* sum_o_STN + ...
        W_CON_GPe_out .* o_CON(:,step) + ...
        W_GPe_out_GPe_out .* o_GPe_out(:,step-1);
    a_GPe_out = (a_GPe_out - u_GPe_out) * decay_constant + u_GPe_out;
    o_GPe_out(:,step) = ramp( a_GPe_out, e_GPe, m );
    
    % GPe - Inner
    u_GPe_inn = W_STN_GPe_inn .* sum_o_STN + ...
        W_CON_GPe_inn .* o_CON(:,step) + ...
        W_GPe_out_GPe_inn .* o_GPe_out(:,step) + ...
        W_GPe_inn_GPe_inn .* o_GPe_inn(:,step-1);
    a_GPe_inn = (a_GPe_inn - u_GPe_inn) * decay_constant + u_GPe_inn;
    o_GPe_inn(:,step) = ramp( a_GPe_inn, e_GPe, m );
    
    % GPe - TA
    u_GPe_ta = W_STN_GPe_ta .* sum_o_STN + ...
        W_CON_GPe_ta .* o_CON(:,step) + ...
        W_GPe_out_GPe_ta .* o_GPe_out(:,step) + ...
        W_GPe_inn_GPe_ta .* o_GPe_inn(:,step) + ...
        W_GPe_ta_GPe_ta .* o_GPe_ta(:,step-1);
    a_GPe_ta = (a_GPe_ta - u_GPe_ta) * decay_constant + u_GPe_ta;
    o_GPe_ta(:,step) = ramp( a_GPe_ta, e_GPe, m );
    
    
    % GPi
    u_GPi = W_STN_GPi .* sum_o_STN +...
        W_SEL_GPi .* o_SEL(:,step) +...
        W_GPe_out_GPi .* o_GPe_out(:,step) +...
        W_GPe_inn_GPi .* o_GPe_inn(:,step);
    a_GPi = (a_GPi - u_GPi) * decay_constant + u_GPi;
    o_GPi(:,step) = ramp( a_GPi, e_GPi, m );
end




