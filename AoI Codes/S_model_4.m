clc;
clear;

% Define the state space components
AoI1_values = (1:5);
AoI2_values = (1:5);
B1_level=(0:5);
B2_level=(0:5);
Ss1=[1 2];
Ss2=[1 2];



% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values)* numel(B1_level)* numel(B2_level) * numel(Ss1)* numel(Ss2);

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
        
                        state_space(idx, :) = [aoi1, aoi2,B1,B2,s1,s2];
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



discount = 0.95;

% Loop through all states to calculate transitions and rewards
nn = zeros(1,num_states);

for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
    next_state11 = state;
    if(next_state11(3)>=1 && next_state11(5)==1)
        next_state11(1) = 1;
        next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
        next_state11(3)=next_state11(3)-1;
 
    else
        next_state11(1) = min(next_state11(1) + 1, max(AoI1_values));
        next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
    end
   
  
   
    % Sensor two transitions and rewards
     next_state21 = state;
     if(next_state21(4)>=1 && next_state21(6)==1)
        next_state21(2) = 1;
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
         next_state21(4)=next_state21(4)-1;
        
     else
     next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
     next_state21(2) = min(next_state21(2) + 1, max(AoI2_values));
     end
      
   
    % Find indices of next states in state_space
    next_idx11 = find(ismember(state_space, next_state11, 'rows'));
    
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));
    
    
    % Update transition probability matrices
    P1(i, next_idx11) = 1;
    
    P2(i, next_idx21) = 1;


   
    
    % Update reward matrices
    R1(i) = 1 / (.9*next_state11(1)+.1*next_state11(2));
    R2(i) = 1 / (.9*next_state21(1)+.1*next_state21(2));
    
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;


% Combine reward matrices for all sensors
R = [R1, R2];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

% Plotting the Policy


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
axis([floor(min(state_space(:, 1))), ceil(max(state_space(:, 1))), floor(min(state_space(:, 2))), ceil(max(state_space(:, 2)))]);
% Create a legend
legend('Action 1', 'Action 2', 'Action 3');

% Display the plot
grid on;


% Initialize variables for simulation
current_state = [1, 1,2,2,1,1]; % Initial state (you can choose any initial state)
num_steps = 100; % Number of simulation steps


% Initialize a variable to store AoI values over time
aoi_history = zeros(num_steps, 6);
aa=[];
% Simulate the system
for step = 1:num_steps
    % Determine the action to take based on the current state and policy
    current_idx = find(ismember(state_space, current_state, 'rows'));
    action = policy(current_idx);
    aa=[aa;action];
    
    % Update AoI values based on the chosen action
    if (action == 1) % Sensor 1 sends
       if(current_state(3)>=1 && current_state(5)==1)
            current_state(1) = 1;
            current_state(2) = min(current_state(2) + 1,max(AoI2_values));
            current_state(3)=current_state(3)-1;
        else
            current_state(1) = min(current_state(1) + 1,max(AoI1_values));
            current_state(2) = min(current_state(2) + 1,max(AoI2_values));
        end
    else % Sensor 2 sends
        if(current_state(4)>=1 && current_state(6)==1)
           current_state(2) = 1;
           current_state(1) = min(current_state(1) + 1,max(AoI1_values));
           current_state(4)=current_state(4)-1;
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
axis([0 10 0 10]);
grid on;


