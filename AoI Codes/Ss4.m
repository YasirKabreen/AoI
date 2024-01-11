clc;
clear;

% Define the state space components
AoI1_values = [1, 2,3,4,5,6,7,8];
AoI2_values = [1, 2,3,4,5,6,7,8];
LocationAGV_values = [1, 2, 3, 4,5];
%Source1State_values = [1,2];
%Source2State_values = [1,2];

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * ...
    numel(LocationAGV_values) * numel(LocationAGV_values)  ;

% Create a matrix to store all combinations of states
state_space = zeros(num_states, 4); % e columns for state components

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for locAGV1 = LocationAGV_values
            for locAGV2 = LocationAGV_values
                
                        state_space(idx, :) = [aoi1, aoi2, locAGV1, locAGV2];
                        idx = idx + 1;
                   
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
    if(state(3)==state(4))||(state(3)+state(4)==3)
        next_state1(2)=1;
    else
    next_state1(2) = min(next_state1(2) + 1, max(AoI2_values));
    end
    nn=[nn ; next_state1];
   
    % Sensor two transitions and rewards
    next_state2 = state;
    next_state2(4) = mod(next_state2(4) + 1, num_zones + 1);
    next_state2(4) = max(next_state2(4), 1); % Ensure it's at least 1
     next_state2(3) = mod(next_state2(3) + 1, num_zones + 1);
    next_state2(3) = max(next_state2(3), 1); % Ensure it's at least 1
     if(state(3)==state(4))||(state(3)+state(4)==3)
        next_state1(1)=1;
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


% Calculate the average AoI
average_AoI = 0;
for i = 1:num_states
    initial_state_prob = 1 / num_states; % Assuming equal probability for each state
    average_AoI = average_AoI + initial_state_prob * V(i);
end

fprintf('Average AoI: %.4f\n', average_AoI);



% Initialize the total AoI
total_AoI = 0;

% Specify the number of time steps for the simulation
num_time_steps = 100; % You can adjust this as needed

% Specify the initial state (choose an appropriate initial state)
initial_state = [1, 1, 1, 1]; % Adjust as needed

% Simulate the system
current_state = initial_state;

for t = 1:num_time_steps
    % Use the policy to determine the action to take
    action = policy(find(ismember(state_space, current_state, 'rows')));
    
    % Transition to the next state based on the chosen action
    if action == 1
        % Perform Sensor 1 action
        % Update current_state accordingly
        % Update AoI for Sensor 1
    elseif action == 2
        % Perform Sensor 2 action
        % Update current_state accordingly
        % Update AoI for Sensor 2
    else
        % No action
        % Update current_state accordingly (if needed)
    end
    
    % Calculate and update the total AoI
    total_AoI = total_AoI + current_state(1) + current_state(2); % Assuming AoI1 and AoI2 are stored in current_state
    
    % Check for termination condition (if needed)
end

fprintf('Total AoI: %d\n', total_AoI);

