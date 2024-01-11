clc;
clear;

% Define the state space components including AoI3
AoI1_values = [1, 2, 3, ];
AoI2_values = [1, 2, 3, ];
AoI3_values = [1, 2, 3, ];
LocationAGV_values = [1, 2, 3, 4, 5];

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(AoI3_values) * numel(LocationAGV_values) * numel(LocationAGV_values)* numel(LocationAGV_values);

% Create a matrix to store all combinations of states with additional AoI3
state_space = zeros(num_states, 6); % Five columns for state components, including AoI3

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for aoi3 = AoI3_values % Added AoI3 loop
            for locAGV1 = LocationAGV_values
                for locAGV2 = LocationAGV_values
                    for locAGV3 = LocationAGV_values
                    state_space(idx, :) = [aoi1, aoi2, aoi3, locAGV1, locAGV2,locAGV3];
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
P3 = zeros(num_states, num_states); % Transition matrix for Sensor 3

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
R3 = zeros(num_states, 1); % Reward matrix for Sensor 3

num_zones = max(LocationAGV_values);
discount = 0.95;

% Loop through all states to calculate transitions and rewards
nn = [];
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    next_state1 = state;
    next_state1(1) = 1;
    next_state1(4) = mod(next_state1(4) + 1, num_zones + 1);
    next_state1(4) = max(next_state1(4), 1); % Ensure it's at least 1
    next_state1(5) = mod(next_state1(5) + 1, num_zones + 1);
    next_state1(5) = max(next_state1(5), 1); % Ensure it's at least 1
    if (state(4) == state(5)) || (state(4) + state(5) == 3)
        next_state1(2) = 1;
    else
        next_state1(2) = min(next_state1(2) + 1, max(AoI2_values));
    end
    nn = [nn; next_state1];
   
    % Sensor two transitions and rewards
    next_state2 = state;
    next_state2(5) = mod(next_state2(5) + 1, num_zones + 1);
    next_state2(5) = max(next_state2(5), 1); % Ensure it's at least 1
    next_state2(4) = mod(next_state2(4) + 1, num_zones + 1);
    next_state2(4) = max(next_state2(4), 1); % Ensure it's at least 1
    if (state(4) == state(5)) || (state(4) + state(5) == 3)
        next_state1(1) = 1;
    else
        next_state2(1) = min(next_state2(1) + 1, max(AoI1_values));
    end
    
    % Sensor three transitions and rewards
    next_state3 = state;
    next_state3(4) = mod(next_state3(4) + 1, num_zones + 1);
    next_state3(4) = max(next_state3(4), 1); % Ensure it's at least 1
    next_state3(5) = mod(next_state3(5) + 1, num_zones + 1);
    next_state3(5) = max(next_state3(5), 1); % Ensure it's at least 1
    if (state(4) == state(5)) || (state(4) + state(5) == 3)
        next_state3(3) = 1;
    else
        next_state3(3) = min(next_state3(3) + 1, max(AoI3_values));
    end
    
    % Find indices of next states in state_space
    next_idx1 = find(ismember(state_space, next_state1, 'rows'));
    next_idx2 = find(ismember(state_space, next_state2, 'rows'));
    next_idx3 = find(ismember(state_space, next_state3, 'rows'));
    
    % Update transition probability matrices
    P1(i, next_idx1) = 1;
    P2(i, next_idx2) = 1;
    P3(i, next_idx3) = 1;
    
    % Update reward matrices
    R1(i) = 1 / state(1);
    R2(i) = 1 / state(2);
    R3(i) = 1 / state(3);
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;

% Combine reward matrices for all sensors
R = [R1, R2, R3];

% Perform policy iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

% Plotting the Policy
figure;
plot(policy);
xlabel('State Index');
ylabel('Action');
title('Optimal Policy for AoI Minimization');
axis([0 length(state_space) + 1 0 4]); % Adjust the axis based on the number of sensors
yticks(0:4);
yticklabels({'No Action', 'Sensor 1', 'Sensor 2', 'Sensor 3'});
grid on;

% Initialize the total AoI
total_AoI = 0;

% Specify the number of time steps for the simulation
num_time_steps = 100; % You can adjust this as needed

% Specify the initial state (choose an appropriate initial state)
initial_state = [1, 1, 1, 4, 2]; % Adjust as needed

% Simulate the system
current_state = initial_state;
vv = [];
for t = 1:num_time_steps
    % Use the policy to determine the action
    action = policy(find(ismember(state_space, current_state, 'rows')));
    
    % Transition to the next state based on the chosen action
    if action == 1
        % Perform Sensor 1 action
        current_state(1) = 1;
        if (current_state(4) == current_state(5)) || (current_state(4) + current_state(5) == 3)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI2_values));
        end
        current_state(3) = min(current_state(3) + 1, max(AoI3_values));
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
    elseif action == 2
        % Perform Sensor 2 action
        current_state(2) = 1;
        if (current_state(4) == current_state(5)) || (current_state(4) + current_state(5) == 3)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        current_state(3) = min(current_state(3) + 1, max(AoI3_values));
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
    elseif action == 3
        % Perform Sensor 3 action
        current_state(3) = 1;
        if (current_state(4) == current_state(5)) || (current_state(4) + current_state(5) == 3)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        current_state(2) = min(current_state(2) + 1, max(AoI2_values));
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
    else
        % No action
        % Example: If no action, you may just update AGV locations
        % Update AGV locations or other state variables
    end
    
    % Calculate and update the total AoI
    total_AoI = total_AoI + current_state(1) + current_state(2) + current_state(3); % Update for all sensors
    
    % Check for termination condition (if needed)
    % For example, you can add a termination condition based on total_AoI or t
    % Example: Terminate if total_AoI exceeds a certain threshold
    if total_AoI > 1000000000
        break; % Terminate the simulation
    end
    vv = [vv; current_state];
end

fprintf('Total AoI: %d\n', total_AoI);


% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);

num_zones = max(LocationAGV_values);
discount = 0.95;

% Loop through all states to calculate transitions and rewards
nn = [];
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    next_state1 = state;
    next_state1(1) = 1;
    next_state1(3) = mod(next_state1(3) + 1, num_zones + 1);
    next_state1(3) = max(next_state1(3), 1); % Ensure it's at least 1
    next_state1(4) = mod(next_state1(4) + 1, num_zones + 1);
    next_state1(4) = max(next_state1(4), 1); % Ensure it's at least 1
    if (state(3) == state(4)) || (state(3) + state(4) == 3)
        next_state1(2) = 1;
    else
        next_state1(2) = min(next_state1(2) + 1, max(AoI2_values));
    end
    nn = [nn; next_state1];
    
    % Sensor two transitions and rewards
    next_state2 = state;
    next_state2(4) = mod(next_state2(4) + 1, num_zones + 1);
    next_state2(4) = max(next_state2(4), 1); % Ensure it's at least 1
    next_state2(3) = mod(next_state2(3) + 1, num_zones + 1);
    next_state2(3) = max(next_state2(3), 1); % Ensure it's at least 1
    if (state(3) == state(4)) || (state(3) + state(4) == 3)
        next_state2(1) = 1;
    else
        next_state2(1) = min(next_state2(1) + 1, max(AoI1_values));
    end
    
    % Find indices of next states in state_space
    next_idx1 = find(ismember(state_space, next_state1, 'rows'));
    next_idx2 = find(ismember(state_space, next_state2, 'rows'));
    
    % Update transition probability matrices
    P1(i, next_idx1) = 1;
    P2(i, next_idx2) = 1;
    
    % Update reward matrices
    R1(i) = 1 / state(1);
    R2(i) = 1 / state(2);
end

% Combine transition probability matrices for both sensors
P(:, :, 1) = P1;
P(:, :, 2) = P2;

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
axis([0 length(state_space) + 1 0 3]);
% xticks(1:length(state_space));
yticks(0:3);
yticklabels({'No Action', 'Sensor 1', 'Sensor 2'});
grid on;

% Initialize the total AoI
total_AoI = 0;

% Specify the number of time steps for the simulation
num_time_steps = 100; % You can adjust this as needed

% Specify the initial state (choose an appropriate initial state)
initial_state = [1,1, 4, 2]; % Adjust as needed

% Simulate the system
current_state = initial_state;
vv = [];
action_sequence = [];

for t = 1:num_time_steps
    % Use the policy to determine the action
    state_idx = find(all(bsxfun(@eq, state_space, current_state), 2));
    action = policy(state_idx);
    
    % Store the chosen action in the action sequence
    action_sequence = [action_sequence; state_idx];
    
    % Transition to the next state based on the chosen action
    if action == 1
        % Perform Sensor 1 action
        current_state(1) = 1;
        if (current_state(3) == current_state(4)) || (current_state(3) + current_state(4) == 3)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI2_values));
    
        end
        current_state(3) = mod(current_state(3) + 1, num_zones + 1);
        current_state(3) = max(current_state(3), 1); % Ensure it's at least 1
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
    elseif action == 2
        % Perform Sensor 2 action
        current_state(2) = 1;
        if (current_state(3) == current_state(4)) || (current_state(3) + current_state(4) == 3)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        % Update AGV locations or other state variables as needed
        current_state(3) = mod(current_state(3) + 1, num_zones + 1);
        current_state(3) = max(current_state(3), 1); % Ensure it's at least 1
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
    else
        % No action
        % Update AGV locations or other state variables as needed
    end
    
    % Calculate and update the total AoI
    total_AoI = total_AoI + current_state(1) + current_state(2); % Assuming AoI1 and AoI2 are stored in current_state
    
   
    vv = [vv; current_state];
end

fprintf('Total AoI: %d\n', total_AoI);
