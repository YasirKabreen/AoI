clc;
clear;

% Define the state space components
AoI1_values = [1, 2,3,4,5,6];
AoI2_values = [1, 2,3,4,5,6];
LocationAGV_values = [1, 2, 3, 4];
Source1State_values = [1];
Source2State_values = [1];

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * ...
    numel(LocationAGV_values) * numel(LocationAGV_values) * ...
    numel(Source1State_values) * numel(Source2State_values);

% Create a matrix to store all combinations of states
state_space = zeros(num_states, 6); % 6 columns for state components

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for locAGV1 = LocationAGV_values
            for locAGV2 = LocationAGV_values
                for source1 = Source1State_values
                    for source2 = Source2State_values
                        state_space(idx, :) = [aoi1, aoi2, locAGV1, locAGV2, source1, source2];
                        idx = idx + 1;
                    end
                end
            end
        end
    end
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);


num_zones = max(LocationAGV_values);
discount = 0.95;

% Loop through all states to calculate transitions and rewards
nn=[];
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    next_state1 = state;
    next_state1(1)=1;
    next_state1(3) = mod(next_state1(3) + 1, num_zones + 1);
    next_state1(3) = max(next_state1(3), 1); % Ensure it's at least 1
    next_state1(4) = mod(next_state1(4) + 1, num_zones + 1);
    next_state1(4) = max(next_state1(4), 1); % Ensure it's at least 1
    next_state1(2) = min(next_state1(2) + 1, max(AoI2_values));
    nn=[nn ; next_state1];
   
    % Sensor two transitions and rewards
    next_state2 = state;
    next_state2(4) = mod(next_state2(4) + 1, num_zones + 1);
    next_state2(4) = max(next_state2(4), 1); % Ensure it's at least 1
     next_state2(3) = mod(next_state2(3) + 1, num_zones + 1);
    next_state2(3) = max(next_state2(3), 1); % Ensure it's at least 1
    next_state2(1) = min(next_state2(1) + 1, max(AoI1_values));
    
    % Find indices of next states in state_space
    next_idx1 = find(ismember(state_space, next_state1, 'rows'));
    next_idx2 = find(ismember(state_space, next_state2, 'rows'));
    
    % Update transition probability matrices
    P1(i, next_idx1) = 1;
    P2(i, next_idx2) = 1;
    
    % Update reward matrices
    R1(i) =  1/state(1);
    R2(i) = 1/ state(2);
end

% Combine transition probability matrices for both sensors
P(:,:,1) = P1;
P(:,:,2) = P2;

% Combine reward matrices for both sensors
R = [R1, R2];

% Perform policy iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

    % Plotting the Policy
figure;
plot(policy);
xlabel('State Index');
ylabel('Action');
title('Optimal Policy for AoI Minimization');
axis([0 length(state_space)+1 0 3]);
%xticks(1:length(state_space));
yticks(0:3);
yticklabels({ 'No Action', 'Sensor 1', 'Sensor 2'});
grid on;
