clc;
clear;

% Define the state space components
AoI1_values = (1:5);
AoI2_values = (1:5);
B1_level = (0:5);
B2_level = (0:5);
Ss1 = [1 2];
Ss2 = [1 2];

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(B1_level) * numel(B2_level) * numel(Ss1) * numel(Ss2);

% Create a matrix to store all combinations of states
state_space = zeros(num_states, 6);

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for B1 = B1_level
            for B2 = B2_level
                for s1 = Ss1
                    for s2 = Ss2
                        state_space(idx, :) = [aoi1, aoi2, B1, B2, s1, s2];
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

% Initialize reward matrices
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);

discount = 0.95;

% Loop through all states to calculate transitions and rewards
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    if state(5) == 1 && state(3) > 0 % Sensor 1 is ON and has battery
        next_state11 = state;
        next_state11(1) = 1; % Reset AoI to 1
        next_state11(3) = max(0, state(3) - 1); % Consume 1 unit of battery
        next_state11(2) = state(2) + 1; % Increase AoI of Sensor 2 by 1
    else
        % Sensor 1 is OFF or has no battery, AoI increases by 1 for both sensors
        next_state11 = state;
        next_state11(1) = state(1) + 1;
        next_state11(2) = state(2) + 1;
    end

    % Sensor two transitions and rewards
    if state(6) == 1 && state(4) > 0 % Sensor 2 is ON and has battery
        next_state21 = state;
        next_state21(2) = 1; % Reset AoI to 1
        next_state21(4) = max(0, state(4) - 1); % Consume 1 unit of battery
        next_state21(1) = state(1) + 1; % Increase AoI of Sensor 1 by 1
    else
        % Sensor 2 is OFF or has no battery, AoI increases by 1 for both sensors
        next_state21 = state;
        next_state21(2) = state(2) + 1;
        next_state21(1) = state(1) + 1;
    end

    % Find indices of next states in state_space
    next_idx11 = find(ismember(state_space, next_state11, 'rows'));
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));

    % Update transition probability matrices
    P1(i, next_idx11) = 1;
    P2(i, next_idx21) = 1;

    % Calculate the AoI difference between the current and next states
    aoi_difference1 = state(1) - next_state11(1);
    aoi_difference2 = state(2) - next_state21(2);

    % Penalize higher AoI values (lower values are better)
    R1(i) = -aoi_difference1;
    R2(i) = -aoi_difference2;
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;

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
yticks(0:4);
yticklabels({'No Action', 'Sensor 1', 'Sensor 2', 'Sensor 3'});
grid on;

% Initialize variables for simulation
current_state = [1, 1, 2, 2, 1, 1]; % Initial state (you can choose any initial state)
num_steps = 100; % Number of simulation steps

% Initialize a variable to store AoI values over time
aoi_history = zeros(num_steps, 6);
aa = [];

% Simulate the system
for step = 1:num_steps
    % Determine the action to take based on the current state and policy
    current_idx = find(ismember(state_space, current_state, 'rows'));
    action = policy(current_idx);
    aa = [aa; action];
    
    % Update AoI values based on the chosen action
    if (action == 1) % Sensor 1 sends
        if (current_state(3) >= 1 && current_state(5) == 1)
            current_state(1) = 1;
            current_state(2) = min(current_state(2) + 1, max(AoI2_values));
            current_state(3) = current_state(3) - 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
            current_state(2) = min(current_state(2) + 1, max(AoI2_values));
        end
    else % Sensor 2 sends
        if (current_state(4) >= 1 && current_state(6) == 1)
            current_state(2) = 1;
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
            current_state(4) = current_state(4) - 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
            current_state(2) = min(current_state(2) + 1, max(AoI2_values));
        end
    end
    
    % Store AoI values in the history
    aoi_history(step, :) = current_state;
end

% Plot the AoI values over time
figure;
plot(1:num_steps, aoi_history(:, 1), 'b-', 1:num_steps, aoi_history(:, 2), 'r-');
xlabel('Time Step');
ylabel('AoI');
legend('Sensor 1 AoI', 'Sensor 2 AoI');
title('AoI Evolution Over Time');
axis([0 10 0 10]);
grid on;
