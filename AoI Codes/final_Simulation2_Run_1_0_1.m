clc;
clear;

% Define the state space components
AoI1_values = (1:30);
B_level1=(0:15);

% Calculate the number of states
num_states = numel(AoI1_values) * numel(B_level1);

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 2); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for b1 = B_level1
        state_space(idx, :) = [aoi1,b1];
        idx = idx + 1;
    end
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);

discount = 0.95;

e= 0.2; % Probability of EH
q= input('enter the probability Source coverage q: '); % Probability of Source coverage
p= 0.9; % Probability of Successful trasnmission

for i = 1:num_states
    state = state_space(i, :);
    
    % Action 1: Sending Fresh Data or Not Sending Data with Empty Update
    
    if state(2) > 0 % If battery is greater than 0
        next_state11 = [1, state(2)-1+1];
        next_state12 = [1, state(2)-1];
        next_state13 = [min(state(1) + 1, max(AoI1_values)), state(2)-1+1];
        next_state14 = [min(state(1) + 1, max(AoI1_values)), state(2)-1];
        
        % Update transition probability matrices
        P1(i, ismember(state_space, next_state11, 'rows')) = q * p *  e   ;              
        P1(i, ismember(state_space, next_state12, 'rows')) = q * p * (1-e)  ;          
        P1(i, ismember(state_space, next_state13, 'rows')) = (( 1- q * p ) *  e)   ;              
        P1(i, ismember(state_space, next_state14, 'rows')) = (( 1- q * p ) * (1-e)) ;        
        
    else % If battery is 0
        % Scenario 1: e (with energy harvesting)
        next_state11 = [min(state(1) + 1, max(AoI1_values)), min(state(2)+1, max(B_level1))];
        
        % Scenario 2: (1-e) (without energy harvesting)
        next_state12 = [min(state(1) + 1, max(AoI1_values)), state(2)];
        
        % Update transition probability matrices for Action 1 with battery = 0
        P1(i, ismember(state_space, next_state11, 'rows')) = e ;
        P1(i, ismember(state_space, next_state12, 'rows')) = 1 - e ;
    end
    
    % Action 2: Not Sending and Conserving Energy
    next_state21 = [min(state(1) + 1, max(AoI1_values)),  min(state(2)+1, max(B_level1)) ];
    next_state22 = [min(state(1) + 1, max(AoI1_values)), state(2)];
    
    % Update transition probability matrices for Action 2
    P2(i, ismember(state_space, next_state21, 'rows')) = e+  P2(i, ismember(state_space, next_state21, 'rows')) ;
    P2(i, ismember(state_space, next_state22, 'rows')) = (1 - e) + P2(i, ismember(state_space, next_state22, 'rows')) ;
    
    % Compute Costs for each state-action pair
    AoI_cost = state(1);
    Energy_cost_A1 = max(B_level1) - state(2);
    Energy_saved_A2 = state(2);
    w1 = 0.9; w2 = 0.1;
    w3 = 0.9; w4 = 0.1;
   if state(2) > 0 % If battery is greater than 0
    Cost_A1 = w1 * (AoI_cost -p*q*state(1))+ w2 * Energy_cost_A1;
    Cost_A2 = w3 * AoI_cost - w4 * Energy_saved_A2;
   else
    Cost_A1 = w1 * (AoI_cost )+ w2 * Energy_cost_A1;
    Cost_A2 = w3 * AoI_cost - w4 * Energy_saved_A2;
   end
    
    % Update reward matrices
    R1(i) = -Cost_A1;
    R2(i) = -Cost_A2;
    
   
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;


% Combine reward matrices for all sensors
R = [R1, R2];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);



% Create a figure
figure;

% Define markers and colors for the policy
markers = {'g^', 'ro','b*'};
colors = {'green', 'red','blue'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(policy)
    action = policy(i);
    marker = markers{action};
    color = colors{action};
    scatter(state_space(i, 2), state_space(i, 1), 50, marker, 'filled', 'MarkerEdgeColor', color);
    hold on; % To ensure points are plotted on the same graph
end

% Add labels
xlabel('Battery level');
ylabel('AoI');
% Set the X and Y axis limits with a step size of 1
axis([0 16 0 31]);
xticks(1:1:35);
yticks(0:1:35);
% Create a legend
legend('Action 1', 'Action 2');

% Display the plot
grid on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters
initial_state = [1, 15]; % You can change this to any initial state of your choice

% Parameters
num_steps = 1000000;

sim_states_opt = zeros(num_steps, 2);
sim_actions_opt = zeros(num_steps, 1);
sim_rewards_opt = zeros(num_steps, 1);
sim_states_opt(1, :) = initial_state;
for step = 1:num_steps
    current_state = sim_states_opt(step, :);
    state_idx = find(ismember(state_space, current_state, 'rows'));
    action = policy(state_idx);
    sim_actions_opt(step) = action;
    next_state_idx = randsample(num_states, 1, true, P(state_idx,:,action));
    next_state = state_space(next_state_idx, :);
    sim_states_opt(step+1, :) = next_state;
    sim_rewards_opt(step) = R(state_idx, action);
end

% Simulate using random policy
sim_states_rand = zeros(num_steps, 2);
sim_actions_rand = zeros(num_steps, 1);
sim_rewards_rand = zeros(num_steps, 1);
sim_states_rand(1, :) = initial_state;
for step = 1:num_steps
    current_state = sim_states_rand(step, :);
    state_idx = find(ismember(state_space, current_state, 'rows'));
    action = randi(2);  % Choose action uniformly at random
    sim_actions_rand(step) = action;
    next_state_idx = randsample(num_states, 1, true, P(state_idx,:,action));
    next_state = state_space(next_state_idx, :);
    sim_states_rand(step+1, :) = next_state;
    sim_rewards_rand(step) = R(state_idx, action);
end

% Simulate using greedy policy
sim_states_greedy = zeros(num_steps, 2);
sim_actions_greedy = zeros(num_steps, 1);
sim_rewards_greedy = zeros(num_steps, 1);
sim_states_greedy(1, :) = initial_state;
greedy_policy = zeros(num_states, 1);
for i = 1:num_states
    [~, greedy_policy(i)] = max(R(i, :));
end
for step = 1:num_steps
    current_state = sim_states_greedy(step, :);
    state_idx = find(ismember(state_space, current_state, 'rows'));
    action = greedy_policy(state_idx);
    sim_actions_greedy(step) = action;
    next_state_idx = randsample(num_states, 1, true, P(state_idx,:,action));
    next_state = state_space(next_state_idx, :);
    sim_states_greedy(step+1, :) = next_state;
    sim_rewards_greedy(step) = R(state_idx, action);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize cost matrices as negative rewards (if not already defined as costs)
C1 = -R1;
C2 = -R2;

% Combine cost matrices for all sensors
C = [C1, C2 ];
% Calculate average running cost over time
avg_cost_opt = -cumsum(sim_rewards_opt)./(1:num_steps)';  % Negative sign to convert rewards to costs
avg_cost_rand = -cumsum(sim_rewards_rand)./(1:num_steps)';  % Negative sign to convert rewards to costs
avg_cost_greedy = -cumsum(sim_rewards_greedy)./(1:num_steps)';  % Negative sign to convert rewards to costs
% Plot the average running cost over time
figure;
plot(avg_cost_opt, 'r', 'LineWidth', 2);
hold on;
plot(avg_cost_rand, 'b', 'LineWidth', 2);
plot(avg_cost_greedy, 'g', 'LineWidth', 2);
xlabel('Time step');
ylabel('Average Cost');
legend('Optimal Policy', 'Random Policy', 'Greedy Policy');
title('Average Cost Evolution Over Time');

