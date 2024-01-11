clc;
clear;

% Define the state space components
AoI1_values = (1:20);
AoI2_values = (1:20);
req1 = [0 1];
req2 = [0 1];

% Define transition probabilities
prob = 0.9;
Re1 = 0.1;
Re2 = 0.9;

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values)* numel(req1)* numel(req2);

% Create a matrix to store all combinations of states
state_space = zeros(num_states, 4);

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for re1 = req1
            for re2 = req2
                state_space(idx, :) = [aoi1, aoi2, re1, re2];
                idx = idx + 1;
            end
        end
    end
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);

% Initialize reward matrices
R1 = 1 ./ (state_space(:, 1) + state_space(:, 2));
R2 = R1; % Same rewards for both sensors

discount = 0.95;

% Loop through all states to calculate transitions and rewards
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    next_state11 = state;
    next_state12 = state;
    
    % Send with probability prob.
    if state(3) == 1 % There is a request
        next_state11(1) = 1;
        next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
        if rand() <= Re1
            next_state11(3) = 1;
        else
            next_state11(3) = 0;
        end
    else
        next_state11(1) = min(next_state11(1) + 1, max(AoI1_values));
        next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
        if rand() <= Re1
            next_state11(3) = 1;
        else
            next_state11(3) = 0;
        end
    end
    % Not send with probability 1-prob.
    next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
    next_state12(2) = min(next_state12(2) + 1, max(AoI2_values));
    if rand() <= Re1
        next_state12(3) = 1;
    else
        next_state12(3) = 0;
    end
    
    % Sensor two transitions and rewards
    next_state21 = state;
    next_state22 = state;
    
    % Send with probability prob.
    if state(4) == 1 % There is a request
        next_state21(2) = 1;
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
        if rand() <= Re2
            next_state21(4) = 1;
        else
            next_state21(4) = 0;
        end
    else
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
        next_state21(2) = min(next_state21(2) + 1, max(AoI2_values));
        if rand() <= Re2
            next_state21(4) = 1;
        else
            next_state21(4) = 0;
        end
    end
    
    % Not send with probability 1-prob.
    next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
    next_state22(2) = min(next_state22(2) + 1, max(AoI2_values));
    if rand() <= Re2
        next_state22(4) = 1;
    else
        next_state22(4) = 0;
    end
    
    % Find indices of next states in state_space
  
next_idx11 = find(ismember(state_space, next_state11, 'rows'));
next_idx12 = find(ismember(state_space, next_state12, 'rows'));
next_idx21 = find(ismember(state_space, next_state21, 'rows'));
next_idx22 = find(ismember(state_space, next_state22, 'rows'));

% Update transition probability matrices
P1(i, next_idx11) = prob;
P1(i, next_idx12) = 1 - prob;
P2(i, next_idx21) = prob;
P2(i, next_idx22) = 1 - prob;

end

% Combine transition probability matrices for all sensors
P = zeros(num_states, num_states, 2);
P(:, :, 1) = P1;
P(:, :, 2) = P2;

% Combine reward matrices for all sensors
R = [R1, R2];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

% Plotting the Policy
figure;
plot(policy);
xlabel('State Index');
ylabel('Action');
title('Optimal Policy for AoI Minimization');
axis([0 length(state_space) + 1 0 4]); % Adjust the axis based on the number of sensors
yticks(0:2);
yticklabels({'No Action', 'Sensor 1', 'Sensor 2'});
grid on;


% Create a figure
figure;

% Define markers and colors for the policy
markers = {'b*', 'ro'}; % Define different markers for Sensor 1 and Sensor 2
colors = {'blue', 'red'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(policy)
    action = policy(i);
    marker = markers{action};
    color = colors{action};
    scatter(state_space(i, 1), state_space(i, 2), 100, marker, 'filled', 'MarkerEdgeColor', color);
    hold on; % To ensure points are plotted on the same graph
end

% Add labels
xlabel('AoI1');
ylabel('AoI2');
% Set the X and Y axis limits with a step size of 1
axis([0 21 0 21]);
xticks(0:1:21);
yticks(0:1:21);
% Create a legend with updated labels
legend('Sensor 1', 'Sensor 2');
