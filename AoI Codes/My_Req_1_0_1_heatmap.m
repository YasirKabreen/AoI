% Define the state space for age of data in cache for both sensors and battery levels
T = 10; % Maximum age of data
B_max = 10; % Maximum battery level
state_space = combvec(1:T, 1:T, 1:B_max, 1:B_max)'; % Cartesian product

% Define possible actions for BS
actions = ["A1", "A2", "A3", "A4"];

% Define user request probabilities for each sensor
p1_S1 = 0.6; 
p2_S1 = 0.1; 
p1_S2 = 0.9; 
p2_S2 = 0.2;

% Initialize reward and transition matrices
R = zeros(size(state_space, 1), length(actions));
P = zeros(size(state_space, 1), size(state_space, 1), length(actions));

% Iterate over state space and compute rewards and transitions
for i = 1:size(state_space, 1)
    current_state = state_space(i, :);
    current_age_S1 = current_state(1);
    current_age_S2 = current_state(2);
    battery_S1 = current_state(3);
    battery_S2 = current_state(4);
    
    for j = 1:length(actions)
        action = actions(j);
        
        switch action
            case "A1"
                % Only age increases
                next_state = min([current_age_S1+1, current_age_S2+1, battery_S1, battery_S2], [T, T, B_max, B_max]);
                P(i, ismember(state_space, next_state, 'rows'), j) = 1;
                R(i, j) = (p1_S1 + p2_S1) / current_age_S1 + (p1_S2 + p2_S2) / current_age_S2;
                
            case "A2"
                % Fresh data from Sensor 1, age of Sensor 2 data increases
                if battery_S1 > 0
                    next_state = [1, min(current_age_S2+1, T), battery_S1-1, battery_S2];
                    P(i, ismember(state_space, next_state, 'rows'), j) = 1;
                    R(i, j) = (p1_S1 + p2_S1);
                else
                    R(i, j) = -inf; % Heavy penalty for attempting action with no battery
                end
                
            case "A3"
                % Fresh data from Sensor 2, age of Sensor 1 data increases
                if battery_S2 > 0
                    next_state = [min(current_age_S1+1, T), 1, battery_S1, battery_S2-1];
                    P(i, ismember(state_space, next_state, 'rows'), j) = 1;
                    R(i, j) = (p1_S2 + p2_S2);
                else
                    R(i, j) = -inf; % Heavy penalty
                end
                
            case "A4"
                % Fresh data from both sensors
                if battery_S1 > 0 && battery_S2 > 0
                    next_state = [1, 1, battery_S1-1, battery_S2-1];
                    P(i, ismember(state_space, next_state, 'rows'), j) = 1;
                    R(i, j) = (p1_S1 + p2_S1) + (p1_S2 + p2_S2);
                else
                    R(i, j) = -inf; % Heavy penalty if any sensor's battery is 0
                end
        end
    end
end

% Compute optimal policy using value iteration
discount = 0.95;
[V, policy] = mdp_policy_iteration(P, R, discount);
% Given values
T = 10;
B = 10; % Replace with your maximum battery level

% Assuming your state_space and policy are already defined

% Choose a specific battery level for visualization
fixed_battery_S1 = 5; % Replace with a specific value
fixed_battery_S2 = 1; % Replace with a specific value

% Filter the policy for the fixed battery levels
filtered_policy = policy(state_space(:,3) == fixed_battery_S1 & state_space(:,4) == fixed_battery_S2);

% Reshape the filtered policy to a matrix form for heatmap
policy_matrix = reshape(filtered_policy, T, T);

% Create a heatmap
figure;
imagesc(policy_matrix);
colorbar;

% Define color map such that 1 -> Red, 2 -> Green, 3 -> Blue, 4 -> Yellow
colormap([1 0 0; 0 1 0; 0 0 1; 1 1 0]);

% Labeling the heatmap
xlabel('Age of Information for Sensor 1');
ylabel('Age of Information for Sensor 2');
title(['Policy Visualization for Battery Levels: ', num2str(fixed_battery_S1), ' & ', num2str(fixed_battery_S2)]);
xticks(1:T);
yticks(1:T);


