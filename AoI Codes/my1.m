clc;
clear;

% Define the state space components
AoI1_values = (1:30);



% Calculate the number of states
num_states = numel(AoI1_values) ;

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 1); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    
                        state_space(idx, :) = aoi1;
                        idx = idx + 1;
       
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);



% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);



discount = 0.95;

% Loop through all states to calculate transitions and rewards

q = 0.6;
for i = 1:num_states
    state = state_space(i, :);
    
       % Sensor one transitions and rewards
    next_state11 = state;
    next_state12 = state;
  
  
        next_state11(1) = 1;
   
   
        next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
       
       

     
   
     % Sensor two transitions and rewards
    next_state21 = state;
    next_state22 = state;
    
        next_state21(1) = 1;
        
        next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
       
     
   
    
   
    % Find indices of next states in state_space
    next_idx11 = find(ismember(state_space, next_state11, 'rows'));
       next_idx12 = find(ismember(state_space, next_state12, 'rows'));
          
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));
       next_idx22 = find(ismember(state_space, next_state22, 'rows'));
        
   
    
    % Update transition probability matrices
    P1(i, next_idx11) = q;
    P1(i, next_idx12) = 1-q;
   
    P2(i, next_idx21) = 1-q;
    P2(i, next_idx22) = q;
   
        
   
   
 

   
    
    % Update reward matrices
    R1(i) = 1 /state(1);
    R2(i) = 1/state(1);
    
    
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

% Create a figure
figure;

% Define markers and colors for the policy
markers = {'b*', 'ro','g^'};
colors = {'blue', 'red','green'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(policy)
    action = policy(i);
    marker = markers{action};
    color = colors{action};
    scatter(state_space(i, 1), state_space(i, 1), 100, marker, 'filled', 'MarkerEdgeColor', color);
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








