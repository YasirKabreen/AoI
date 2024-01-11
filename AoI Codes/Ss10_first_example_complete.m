clc;
clear;

% Define the state space components
AoI1_values = (1:20);
AoI2_values = (1:20);



% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) ;

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

% Loop through all states to calculate transitions and rewards
nn = zeros(1,num_states);
p = input('enter the probability p'); % Probability of sending
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    next_state11 = state;
    next_state12 = state;
    
    %send with probabilty p.
        next_state11(1) = 1;
        next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
     %not send with probabilty 1-p.
     next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
     next_state12(2) = min(next_state12(2) + 1, max(AoI2_values));
    
   
   
   
    % Sensor two transitions and rewards
  next_state21 = state;
    next_state22 = state;
    
    %send with probabilty p.
        next_state21(2) = 1;
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
     %not send with probabilty 1-p.
     next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
     next_state22(2) = min(next_state22(2) + 1, max(AoI2_values));
     
     % Sensor three transitions and rewards
    next_state31 = state;
    next_state32 = state;
    
    %send with probabilty 1-p.
        next_state31(1) = 1;
        next_state31(2) =1;
     %not send with probabilty 1-p.
     next_state32(1) = min(next_state32(1) + 1, max(AoI1_values));
     next_state32(2) = min(next_state32(2) + 1, max(AoI2_values));
    
   
    % Find indices of next states in state_space
    next_idx11 = find(ismember(state_space, next_state11, 'rows'));
     next_idx12 = find(ismember(state_space, next_state12, 'rows'));
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));
    next_idx22 = find(ismember(state_space, next_state22, 'rows'));
     next_idx31 = find(ismember(state_space, next_state31, 'rows'));
      next_idx32 = find(ismember(state_space, next_state32, 'rows'));
    
    % Update transition probability matrices
    P1(i, next_idx11) = p;
    P1(i, next_idx12) = 1-p;
    P2(i, next_idx21) = p;
    P2(i, next_idx22) = 1-p;
    P3(i, next_idx31) = 1-p;
    P3(i, next_idx32) = p;
 

   
    
    % Update reward matrices
    R1(i) = 1 / (state(1)+state(2));
    R2(i) = 1 / (state(1)+state(2));
    R3(i) = 1 / (state(1)+state(2));
    
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;


% Combine reward matrices for all sensors
R = [R1, R2,R3];

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

% Create a figure
figure;

% Define markers and colors for the policy
markers = {'b*', 'ro'};
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
% Create a legend
legend('Action 1', 'Action 2');

% Display the plot
grid on;




% Initialize variables for simulation
current_state = [1, 1]; % Initial state (you can choose any initial state)
num_steps = 100; % Number of simulation steps


% Initialize a variable to store AoI values over time
aoi_history = zeros(num_steps, 2);

% Simulate the system
for step = 1:num_steps
    % Determine the action to take based on the current state and policy
    current_idx = find(ismember(state_space, current_state, 'rows'));
    action = policy(current_idx);
    
    % Update AoI values based on the chosen action
    if action == 1 % Sensor 1 sends
        if rand() <= p % Send with probability p
            current_state(1) = 1;
            current_state(2) = min(current_state(2) + 1,max(AoI2_values));
        else
            current_state(1) = min(current_state(1) + 1,max(AoI1_values));
             current_state(2) = min(current_state(2) + 1,max(AoI2_values));
        end
    elseif action == 2 % Sensor 2 sends
        if rand() <= p % Send with probability p
            current_state(2) = 1;
             current_state(1) = min(current_state(1) + 1,max(AoI1_values));
        else
          current_state(1) = min(current_state(1) + 1,max(AoI1_values));
           current_state(2) = min(current_state(2) + 1,max(AoI2_values));
        end
        elseif action == 3 % Sensor 3 sends
        if rand() <= 1-p % Send with probability 1-p
            current_state(2) = 1;
             current_state(1) = 1;
        else
          current_state(1) = min(current_state(1) + 1,max(AoI1_values));
           current_state(2) = min(current_state(2) + 1,max(AoI2_values));
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
axis([0 100 0 10]);
grid on;

% Assuming you have already defined state_space and policy

% Create a figure for the heatmaps
figure;

% Define the AoI values as cell arrays for labels
AoI1_labels = cellstr(num2str(AoI1_values'));
AoI2_labels = cellstr(num2str(AoI2_values'));

% Create a heatmap for Sensor 1's action
subplot(1, 2, 1);
heatmap(AoI1_labels, AoI2_labels, reshape(policy(:, 1), numel(AoI1_values), numel(AoI2_values)));
title('Sensor 1 Action Policy');
xlabel('AoI1');
ylabel('AoI2');

% Create a heatmap for Sensor 2's action
subplot(1, 2, 2);
heatmap(AoI1_labels, AoI2_labels, reshape(policy(:, 2), numel(AoI1_values), numel(AoI2_values)));
title('Sensor 2 Action Policy');
xlabel('AoI1');
ylabel('AoI2');

% Adjust subplot layout
sgtitle('MDP Policy Heatmaps for Two Sensors');



