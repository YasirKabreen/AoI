clc;
clear;

% Define the state space components
AoI1_values = (1:30);
AoI2_values = (1:30);

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values);

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 2); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        state_space(idx, :) = [aoi1, aoi2];
        idx = idx + 1;
    end
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);
P3 = zeros(num_states, num_states);

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
R3 = zeros(num_states, 1);

discount = 0.95;
q = 1;
p = input('Enter the probability p: '); % Probability of sending

for i = 1:num_states
    state = state_space(i, :);
    delta1 = state(1);
    delta2 = state(2);
    
    % Cost is the sum of AoI values
    cost = delta1 + delta2;
    
    % Rewards (negative of cost)
    R1(i) = -cost;
    R2(i) = -cost;
    R3(i) = -cost;

    % New states for successful transmissions
    new_state1 = [1, min(delta2+1, 20)];  % After action 1
    new_state2 = [min(delta1+1, 20), 1];  % After action 2
    new_state3 = [1, 1];                  % After action 3
    
    % New state for unsuccessful transmission
    new_state_fail = [min(delta1+1, 20), min(delta2+1, 20)]; 

    % Sensor one transitions
    idx1 = find(ismember(state_space, new_state1, 'rows'));
    P1(i, idx1) = p*q;  % Transition probability given action 1 is taken
    
    idx_fail1 = find(ismember(state_space, new_state_fail, 'rows'));
    P1(i, idx_fail1) = 1 - p*q;  % Transition probability given action 1 fails

    % Sensor two transitions
    idx2 = find(ismember(state_space, new_state2, 'rows'));
    P2(i, idx2) = p*q;  % Transition probability given action 2 is taken
    
    idx_fail2 = find(ismember(state_space, new_state_fail, 'rows'));
    P2(i, idx_fail2) = 1 - p*q;  % Transition probability given action 2 fails
    
    % Sensor three transitions
    idx3 = find(ismember(state_space, new_state3, 'rows'));
    P3(i, idx3) = (1-p)*q;  % Transition probability given action 3 is taken
    
    idx_fail3 = find(ismember(state_space, new_state_fail, 'rows'));
    P3(i, idx_fail3) = 1 - (1-p)*q;  % Transition probability given action 3 fails
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;

% Combine reward matrices for all sensors
R = [R1, R2, R3];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

% Visualization and plotting remain the same as in the previous snippet...
% Create a figure
figure;

% Define markers and colors for the policy
markers = {'ro','b*', 'g^'};
colors = {'red','blue', 'green'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(policy)
    action = policy(i);
    marker = markers{action};
    color = colors{action};
    scatter(state_space(i, 1), state_space(i, 2), 50, marker, 'filled', 'MarkerEdgeColor', color);
    hold on;
end

% Add labels
xlabel('AoI1');
ylabel('AoI2');
axis([0 21 0 21]);
xticks(0:1:21);
yticks(0:1:21);
legend('Action 1', 'Action 2', 'Action 3');
grid on;

