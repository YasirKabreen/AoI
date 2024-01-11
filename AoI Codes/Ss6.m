clc;
clear;

% Define the state space components
AoI1_values = [1, 2, 3,4,5];
AoI2_values = [1, 2, 3,4,5];
AoI3_values = [1, 2, 3,4,5];
LocationAGV_values = [1, 2, 3, 4, 5];

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(AoI3_values) * numel(LocationAGV_values) * numel(LocationAGV_values) * numel(LocationAGV_values);

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 6); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for aoi3 = AoI3_values 
            for locAGV1 = LocationAGV_values
                for locAGV2 = LocationAGV_values
                    for locAGV3 = LocationAGV_values
                        state_space(idx, :) = [aoi1, aoi2, aoi3, locAGV1, locAGV2, locAGV3];
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
P3 = zeros(num_states, num_states); 

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
R3 = zeros(num_states, 1); 

num_zones = max(LocationAGV_values);
discount = 0.95;

% Loop through all states to calculate transitions and rewards
nn = [];
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    next_state1 = state;
    if (next_state1(4) == 1 || next_state1(4) == 2)
        next_state1(1) = 1;
    else
        next_state1(1) = min(next_state1(1) + 1, max(AoI1_values));
    end
    if (next_state1(5) == 1 || next_state1(5) == 2)
        next_state1(2) = 1;
    else
        next_state1(2) = min(next_state1(2) + 1, max(AoI1_values));
    end
    if (next_state1(6) == 1 || next_state1(6) == 2)
        next_state1(3) = 1;
    else
        next_state1(3) = min(next_state1(3) + 1, max(AoI1_values));
    end
    
    next_state1(4) = mod(next_state1(4) + 1, num_zones + 1);
    next_state1(4) = max(next_state1(4), 1); % Ensure it's at least 1
    next_state1(5) = mod(next_state1(5) + 1, num_zones + 1);
    next_state1(5) = max(next_state1(5), 1); % Ensure it's at least 1
    next_state1(6) = mod(next_state1(6) + 1, num_zones + 1);
    next_state1(6) = max(next_state1(6), 1); % Ensure it's at least 1
   
    nn = [nn; next_state1];
   
    % Sensor two transitions and rewards
    next_state2 = state;
    if (next_state2(4) == 2)
        next_state2(1) = 1;
    else
        next_state2(1) = min(next_state2(1) + 1, max(AoI1_values));
    end
    if (next_state2(5) == 2)
        next_state2(2) = 1;
    else
        next_state2(2) = min(next_state2(2) + 1, max(AoI1_values));
    end
    if (next_state2(6) == 2)
        next_state2(3) = 1;
    else
        next_state2(3) = min(next_state2(3) + 1, max(AoI1_values));
    end
    
    next_state2(4) = mod(next_state2(4) + 1, num_zones + 1);
    next_state2(4) = max(next_state2(4), 1); % Ensure it's at least 1
    next_state2(5) = mod(next_state2(5) + 1, num_zones + 1);
    next_state2(5) = max(next_state2(5), 1); % Ensure it's at least 1
    next_state2(6) = mod(next_state2(6) + 1, num_zones + 1);
    next_state2(6) = max(next_state2(6), 1); % Ensure it's at least 1
    
    % Sensor three transitions and rewards
    next_state3 = state;
    if (next_state3(4) == 4)
        next_state3(1) = 1;
    else
        next_state3(1) = min(next_state3(1) + 1, max(AoI1_values));
    end
    if (next_state3(5) == 4)
        next_state3(2) = 1;
    else
        next_state3(2) = min(next_state3(2) + 1, max(AoI1_values));
    end
    if (next_state3(6) == 4)
        next_state3(3) = 1;
    else
        next_state3(3) = min(next_state3(3) + 1, max(AoI1_values));
    end
    
    next_state3(4) = mod(next_state3(4) + 1, num_zones + 1);
    next_state3(4) = max(next_state3(4), 1); % Ensure it's at least 1
    next_state3(5) = mod(next_state3(5) + 1, num_zones + 1);
    next_state3(5) = max(next_state3(5), 1); % Ensure it's at least 1
    next_state3(6) = mod(next_state3(6) + 1, num_zones + 1);
    next_state3(6) = max(next_state3(6), 1); % Ensure it's at least 1
    
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

% Perform value iteration to find the optimal value function and policy
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
num_time_steps = 100; 

% Specify the initial state 
initial_state = [1, 1, 4, 1, 2, 3]; 

% Simulate the system
current_state = initial_state;
vv = [];
for t = 1:num_time_steps
    % Use the policy to determine the action to take
    action = policy(find(all(bsxfun(@eq, state_space, current_state), 2), 1));

    
    % Transition to the next state based on the chosen action
    if action == 1
        % Perform Sensor 1 action
        if (current_state(4) == 1 || current_state(4) == 2)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        if (current_state(5) == 1 || current_state(5) == 2)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI1_values));
        end
        if (current_state(6) == 1 || current_state(6) == 2)
            current_state(3) = 1;
        else
            current_state(3) = min(current_state(3) + 1, max(AoI1_values));
        end
        
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
        current_state(6) = mod(current_state(6) + 1, num_zones + 1);
        current_state(6) = max(current_state(6), 1); % Ensure it's at least 1
    elseif action == 2
        % Perform Sensor 2 action
        if (current_state(4) == 2)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        if (current_state(5) == 2)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI1_values));
        end
        if (current_state(6) == 2)
            current_state(3) = 1;
        else
            current_state(3) = min(current_state(3) + 1, max(AoI1_values));
        end
        
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
        current_state(6) = mod(current_state(6) + 1, num_zones + 1);
        current_state(6) = max(current_state(6), 1); % Ensure it's at least 1
    elseif action == 3
        % Perform Sensor 3 action
        if (current_state(4) == 4)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        if (current_state(5) == 4)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI1_values));
        end
        if (current_state(6) == 4)
            current_state(3) = 1;
        else
            current_state(3) = min(current_state(3) + 1, max(AoI1_values));
        end
        
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
        current_state(6) = mod(current_state(6) + 1, num_zones + 1);
        current_state(6) = max(current_state(6), 1); % Ensure it's at least 1
    end
    
    % Calculate and update the total AoI
    total_AoI = total_AoI + current_state(1) + current_state(2) + current_state(3); % Assuming AoI1, AoI2, and AoI3 are stored in current_state
    
   
    
    vv = [vv; current_state];
end

fprintf('Total AoI: %d\n', total_AoI);
