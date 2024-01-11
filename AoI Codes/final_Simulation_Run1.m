clc;
clear;

% Define the state space components
AoI1_values = 1:3;
AoI2_values = 1:3;
B = 3; % Maximum battery level
Battery_values = 0:B;

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(Battery_values)^3;

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 5); 

% Generate all combinations of states
idx = 1;
for b1 = Battery_values
    for b2 = Battery_values
        for b3 = Battery_values
            for aoi1 = AoI1_values
                for aoi2 = AoI2_values
                    state_space(idx, :) = [aoi1, aoi2, b1, b2, b3];
                    idx = idx + 1;
                end
            end
        end
    end
end

% Initialize transition probability matrices and reward matrix
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);
P3 = zeros(num_states, num_states);
P4 = zeros(num_states, num_states);

R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
R3 = zeros(num_states, 1);
R4 = zeros(num_states, 1);
discount=0.95;
p = input('Enter the probability p: ');
q = 1;
p_eh = 0.2; % Probability of energy harvesting
w1=0.9;
w2=0.1;
empty_penalty = 1e3; % High penalty for empty battery

for i = 1:num_states
    state = state_space(i, :);
    delta1 = state(1);
    delta2 = state(2);
    b1 = state(3);
    b2 = state(4);
    b3 = state(5);

    % Battery penalty per sensor
if b1 == 0
    battery_penalty1 = empty_penalty;
else
    battery_penalty1 = b1 * 0.05;
end

if b2 == 0
    battery_penalty2 = empty_penalty;
else
    battery_penalty2 = b2 * 0.05;
end

if b3 == 0
    battery_penalty3 = empty_penalty;
else
    battery_penalty3 = b3 * 0.05;
end
Energy_cost1= B-b1;
Energy_cost2= B-b2;
Energy_cost3= B-b3;
Energy_saved=(b1+b2+b3);

    % Cost computation with battery penalty
    Cost1 = w1*((delta1 + delta2)/2 + 1 - 0.5 * p * q * delta1) + w2*(battery_penalty1+Energy_cost1);
    Cost2 = w1*((delta1 + delta2)/2 + 1 - 0.5 * p * q * delta2) + w2*(battery_penalty2+Energy_cost2);
    Cost3 = w1*((delta1 + delta2)/2 + 1 - 0.5 * (1-p) * q * (delta1 + delta2)) + w2*(battery_penalty3+Energy_cost3);
    Cost4 = w1 * ((delta1 + delta2)/2 + 1) - w2 * Energy_saved ;
    % Rewards
    R1(i) = -Cost1;
    R2(i) = -Cost2;
    R3(i) = -Cost3;
    R4(i) = -Cost4;

   % New states after considering energy harvesting for each action:

% Define the maximum AoI values for clarity
max_AoI1 = max(AoI1_values);
max_AoI2 = max(AoI2_values);

% For successful transmissions with energy harvesting:
new_state1_eh = [1, min(delta2+1, max_AoI2), min(max(b1-1, 0)+1, B), min(b2+1, B), min(b3+1, B)];
new_state2_eh = [min(delta1+1, max_AoI1), 1, min(b1+1, B), min(max(b2-1, 0)+1, B), min(b3+1, B)];
new_state3_eh = [1, 1, min(b1+1, B), min(b2+1, B), min(max(b3-1, 0)+1, B)];

% For unsuccessful transmissions with energy harvesting:
new_state_fail1_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), min(max(b1-1, 0)+1, B), min(b2+1, B), min(b3+1, B)];
new_state_fail2_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), min(b1+1, B), min(max(b2-1, 0)+1, B), min(b3+1, B)];
new_state_fail3_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), min(b1+1, B), min(b2+1, B), min(max(b3-1, 0)+1, B)];

% For successful transmissions without energy harvesting:
new_state1_no_eh = [1, min(delta2+1, max_AoI2), max(b1-1, 0), b2, b3];
new_state2_no_eh = [min(delta1+1, max_AoI1), 1, b1, max(b2-1, 0), b3];
new_state3_no_eh = [1, 1, b1, b2, max(b3-1, 0)];

% For unsuccessful transmissions without energy harvesting:
new_state_fail1_no_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), max(b1-1, 0), b2, b3];
new_state_fail2_no_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), b1, max(b2-1, 0), b3];
new_state_fail3_no_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), b1, b2, max(b3-1, 0)];

% For not sending with energy harvesting:
new_state_not_sending_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), min(b1+1, B), min(b2+1, B), min(b3+1, B)];
% For not sending without energy harvesting:
new_state_not_sending_no_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), b1, b2, b3];

% Finding the index for each new state in the state space after EH:
idx1_eh = find(ismember(state_space, new_state1_eh, 'rows'));
idx2_eh = find(ismember(state_space, new_state2_eh, 'rows'));
idx3_eh = find(ismember(state_space, new_state3_eh, 'rows'));

idx_fail1_eh = find(ismember(state_space, new_state_fail1_eh, 'rows'));
idx_fail2_eh = find(ismember(state_space, new_state_fail2_eh, 'rows'));
idx_fail3_eh = find(ismember(state_space, new_state_fail3_eh, 'rows'));

% Finding the index for each new state in the state space without EH:
idx1_no_eh = find(ismember(state_space, new_state1_no_eh, 'rows'));
idx2_no_eh = find(ismember(state_space, new_state2_no_eh, 'rows'));
idx3_no_eh = find(ismember(state_space, new_state3_no_eh, 'rows'));

idx_fail1_no_eh = find(ismember(state_space, new_state_fail1_no_eh, 'rows'));
idx_fail2_no_eh = find(ismember(state_space, new_state_fail2_no_eh, 'rows'));
idx_fail3_no_eh = find(ismember(state_space, new_state_fail3_no_eh, 'rows'));

% Not sending Finding the index for each new state in the state space :
idx_not_sending_eh = find(ismember(state_space, new_state_not_sending_eh, 'rows'));
idx_not_sending_no_eh = find(ismember(state_space, new_state_not_sending_no_eh, 'rows'));

% Updating the transition probabilities for each state considering EH:
P1(i, idx1_eh) = p * q * p_eh;
P1(i, idx_fail1_eh) = (1 - p*q) * p_eh;

P2(i, idx2_eh) = p * q * p_eh;
P2(i, idx_fail2_eh) = (1 - p*q) * p_eh;

P3(i, idx3_eh) = (1-p) * q * p_eh;
P3(i, idx_fail3_eh) = (1 - (1-p)*q) * p_eh;

% Updating the transition probabilities for each state considering no EH:
P1(i, idx1_no_eh) = p * q * (1 - p_eh);
P1(i, idx_fail1_no_eh) = (1 - p*q) * (1 - p_eh);

P2(i, idx2_no_eh) = p * q * (1 - p_eh);
P2(i, idx_fail2_no_eh) = (1 - p*q) * (1 - p_eh);

P3(i, idx3_no_eh) = (1-p) * q * (1 - p_eh);
P3(i, idx_fail3_no_eh) = (1 - (1-p)*q) * (1 - p_eh);

% Updating the transition probabilities for each state considering no EH:
P4(i, idx_not_sending_eh) =    P4(i, idx_not_sending_eh)+p_eh;
P4(i, idx_not_sending_no_eh) = P4(i, idx_not_sending_no_eh)+ 1 - p_eh;

end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;
P(:,:,4) = P4;

% Combine reward matrices for all sensors
R = [R1, R2, R3,R4];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

% Simulation parameters
initial_state = [3, 3, 2, 2, 2]; % You can change this to any initial state of your choice
num_steps = 1000;

% Initialize matrices to store simulation results
sim_states = zeros(num_steps, 5);
sim_actions = zeros(num_steps, 1);
sim_rewards = zeros(num_steps, 1);

% Set the initial state for the simulation
sim_states(1, :) = initial_state;

for step = 1:num_steps
    current_state = sim_states(step, :);
    
    % Find the index of the current state in the state space
    state_idx = find(ismember(state_space, current_state, 'rows'));
    
    % Use the policy to determine the action
    action = policy(state_idx);
    sim_actions(step) = action;

    % Use the transition probabilities to find the next state
    % Randomly draw the next state based on the transition probabilities
    next_state_idx = randsample(num_states, 1, true, P(state_idx,:,action));
    next_state = state_space(next_state_idx, :);

    sim_states(step+1, :) = next_state;
    sim_rewards(step) = R(state_idx, action);
end

% You can now analyze the results in sim_states, sim_actions, and sim_rewards
% For instance, to get the total reward:
total_reward = sum(sim_rewards);

% To plot the AoI values over time for the two sensors:
figure;
plot(sim_states(:,1));
hold on;
plot(sim_states(:,2));
xlabel('Time step');
ylabel('AoI');
legend('Sensor 1', 'Sensor 2');
axis([0 100 0 4]);
title('AoI evolution over time');

% Similarly, you can plot the battery levels, actions taken, etc.


