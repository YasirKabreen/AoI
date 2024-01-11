clc;
clear;

% Define the state space components
AoI_values = 1:30;
B = 15; % Maximum battery level
max_AoI = max(AoI_values); % Define the maximum AoI value for clarity
Battery_values = 0:B;

% Calculate the number of states
num_states = numel(AoI_values) * numel(Battery_values);

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 2); 

% Generate all combinations of states
idx = 1;
for b = Battery_values
    for aoi = AoI_values
        state_space(idx, :) = [aoi, b];
        idx = idx + 1;
    end
end

% Initialize transition probability matrices and reward matrix
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);

discount = 0.95;
p = input('Enter the probability p: ');
q = 1;
p_eh = input('Enter the probability p_EH: '); % Probability of energy harvesting
    w1 = 0.9; w2 = 0.1;
    w3 = 0.9; w4 = 0.1;
empty_penalty = 1e3; % High penalty for empty battery

for i = 1:num_states
    state = state_space(i, :);
    delta = state(1);
    b = state(2);

    % Battery penalty
    if b == 0
        battery_penalty = empty_penalty;
    else
        battery_penalty =0.05*b ;
    end

    % Calculate the energy cost for action A1 and energy saved for action A2
    Energy_cost_A1 = B - b;
    Energy_saved_A2 = b;

    % AoI cost for Action A1 (sending)
    AoI_cost_A1 = delta + 1 - 0.5 * p * q * delta; 

    % AoI cost for Action A2 (not sending)
    AoI_cost_A2 = delta + 1;

    % Cost computation using the weighted approach (added battery penalty)
    Cost_A1 = w1 * AoI_cost_A1 + w2 *  Energy_cost_A1+battery_penalty;
    Cost_A2 = w3 * AoI_cost_A2 - w4 * Energy_saved_A2 ;

    % Rewards
    R1(i) = -Cost_A1;
    R2(i) = -Cost_A2;
    
    
 

    % For successful transmissions with energy harvesting:
    new_state_eh = [1, min(max(b-1, 0)+1, B)];
    % For unsuccessful transmissions with energy harvesting:
    new_state_fail_eh = [min(delta+1, max_AoI), min(max(b-1, 0)+1, B)];
    % For successful transmissions without energy harvesting:
    new_state_no_eh = [1, max(b-1, 0)];
    % For unsuccessful transmissions without energy harvesting:
    new_state_fail_no_eh = [min(delta+1, max_AoI), max(b-1, 0)];

    
    % For not sending with energy harvesting:
    new_state_not_sending_eh = [min(delta+1, max_AoI), min(b+1, B)];
   
    % For not sending without energy harvesting:
    new_state_not_sending_no_eh = [min(delta+1, max_AoI), b];
    
    % Finding the index for each new state in the state space after EH:
    idx_eh = find(ismember(state_space, new_state_eh, 'rows'));
    idx_fail_eh = find(ismember(state_space, new_state_fail_eh, 'rows'));

    % Finding the index for each new state in the state space without EH:
    idx_no_eh = find(ismember(state_space, new_state_no_eh, 'rows'));
    idx_fail_no_eh = find(ismember(state_space, new_state_fail_no_eh, 'rows'));

    % Not sending Finding the index for each new state in the state space :
    idx_not_sending_eh = find(ismember(state_space, new_state_not_sending_eh, 'rows'));
    idx_not_sending_no_eh = find(ismember(state_space, new_state_not_sending_no_eh, 'rows'));
    
    % Updating the transition probabilities for each state considering EH:
    P1(i, idx_eh) = p * q * p_eh;
    P1(i, idx_fail_eh) = (1 - p*q) * p_eh;

    % Updating the transition probabilities for each state considering no EH:
    P1(i, idx_no_eh) = p * q * (1 - p_eh);
    P1(i, idx_fail_no_eh) = (1 - p*q) * (1 - p_eh);
    
    % Updating the transition probabilities for each state considering no EH:
    P2(i, idx_not_sending_eh) =    P2(i, idx_not_sending_eh)+p_eh;
    P2(i, idx_not_sending_no_eh) = P2(i, idx_not_sending_no_eh)+ 1 - p_eh;
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


