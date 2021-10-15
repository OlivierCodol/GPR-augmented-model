function [c, o_GPi] = dynamics_sim(W_D1, W_D2, Rw)

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
a_GPe = zeros( n_channel , 1 );
a_GPi = zeros( n_channel , 1 );

% output arrays
o_SEL = zeros( n_channel , n_steps );
o_CON = zeros( n_channel , n_steps );
o_STN = zeros( n_channel , n_steps );
o_GPe = zeros( n_channel , n_steps );
o_GPi = zeros( n_channel , n_steps );

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
    u_SEL = c(:, step) * D1_scaler;
    a_SEL = (a_SEL - u_SEL) * decay_constant + u_SEL;
    o_SEL(:,step) = ramp(a_SEL,e_SEL,m);
    
    % STRIATUM D2
    u_CONT = c(:, step) * D2_scaler;
    a_CON = (a_CON - u_CONT) * decay_constant + u_CONT;
    o_CON(:,step) = ramp(a_CON,e_CON,m);
    
    % STN
    u_STN = c(:, step)  * W_STN + W_GPe_STN * o_GPe(:,step-1);
    a_STN = (a_STN - u_STN) * decay_constant + u_STN;
    o_STN(:,step) = ramp(a_STN,e_STN,m);
    
    sum_o_STN = sum(o_STN(:,step));
    
    % GPe
    u_GPe =  sum_o_STN * W_STN_GPe + W_CON_GPe * o_CON(:,step);
    a_GPe = (a_GPe - u_GPe) * decay_constant + u_GPe;
    o_GPe(:,step) = ramp(a_GPe,e_GPe,m);
    
    % GPi
    u_GPi = sum_o_STN * W_STN_GPi + W_GPe_GPi * o_GPe(:,step) + W_SEL_GPi * o_SEL(:,step);
    a_GPi = (a_GPi - u_GPi) * decay_constant + u_GPi;
    o_GPi(:,step) = ramp(a_GPi,e_GPi,m);
end




