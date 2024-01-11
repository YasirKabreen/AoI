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

p = input('enter the probability p'); % Probability of sending
q = input('enter the probability q'); % Probability of source state
for i = 1:num_states
    state = state_space(i, :);
    
       % Sensor one transitions and rewards
    next_state11 = state;
    next_state12 = state;
    next_state13 = state;
    next_state14 = state;
    
    %send with probabilty p and source state q.
        next_state11(1) = 1;
        next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
        
    %not send with probabilty 1-p and source state q.
        next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
        next_state12(2) = min(next_state12(2) + 1, max(AoI2_values));
     
    %send with probabilty p and source state 1-q.
       next_state13(1) = min(next_state13(1) + 1, max(AoI1_values));
       next_state13(2) = min(next_state13(2) + 1, max(AoI2_values));
       
    %not send with probabilty 1-p and source state 1-q.
      next_state14(1) = min(next_state14(1) + 1, max(AoI1_values));
      next_state14(2) = min(next_state14(2) + 1, max(AoI2_values));
   
     % Sensor two transitions and rewards
    next_state21 = state;
    next_state22 = state;
    next_state23 = state;
    next_state24 = state;
    
    %send with probabilty p and source state q.
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
        next_state21(2) = min(next_state21(2) + 1, max(AoI2_values));
        
    %not send with probabilty 1-p and source state q.
        next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
        next_state22(2) = min(next_state22(2) + 1, max(AoI2_values));
     
    %send with probabilty p and source state 1-q.
       next_state23(1) = min(next_state23(1) + 1, max(AoI1_values));
       next_state23(2) = 1;
       
    %not send with probabilty 1-p and source state 1-q.
      next_state24(1) = min(next_state24(1) + 1, max(AoI1_values));
      next_state24(2) = min(next_state24(2) + 1, max(AoI2_values));
     
     % Sensor three transitions and rewards
    next_state31 = state;
    next_state32 = state;
    next_state33 = state;
    next_state34 = state;
    
    %send with probabilty p and source state q.
        next_state31(1) = min(next_state31(1) + 1, max(AoI1_values));
        next_state31(2) = min(next_state31(2) + 1, max(AoI2_values));
        
    %not send with probabilty 1-p and source state q.
        next_state32(1) = 1;
        next_state32(2) = min(next_state32(2) + 1, max(AoI2_values));
     
    %send with probabilty p and source state 1-q.
       next_state33(1) = min(next_state33(1) + 1, max(AoI1_values));
       next_state33(2) = min(next_state33(2) + 1, max(AoI2_values));
       
    %not send with probabilty 1-p and source state 1-q.
      next_state34(1) = min(next_state34(1) + 1, max(AoI1_values));
      next_state34(2) = 1;
    
   
    % Find indices of next states in state_space
    next_idx11 = find(ismember(state_space, next_state11, 'rows'));
       next_idx12 = find(ismember(state_space, next_state12, 'rows'));
          next_idx13 = find(ismember(state_space, next_state13, 'rows'));
             next_idx14 = find(ismember(state_space, next_state14, 'rows'));
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));
       next_idx22 = find(ismember(state_space, next_state22, 'rows'));
          next_idx23 = find(ismember(state_space, next_state23, 'rows'));
             next_idx24 = find(ismember(state_space, next_state24, 'rows'));
    next_idx31 = find(ismember(state_space, next_state31, 'rows'));
       next_idx32 = find(ismember(state_space, next_state32, 'rows'));
          next_idx33 = find(ismember(state_space, next_state33, 'rows'));
             next_idx34 = find(ismember(state_space, next_state34, 'rows'));
    
    % Update transition probability matrices
    P1(i, next_idx11) = p*q;
    P1(i, next_idx12) = (1-p)*q;
    P1(i, next_idx13) = p*(1-q);
    P1(i, next_idx14) = (1-p)*(1-q);
    
    P2(i, next_idx21) = p*q;
    P2(i, next_idx22) = (1-p)*q;
    P2(i, next_idx23) = p*(1-q);
    P2(i, next_idx24) = (1-p)*(1-q);
        
    P3(i, next_idx31) = p*q;
    P3(i, next_idx32) = (1-p)*q;
    P3(i, next_idx33) = p*(1-q);
    P3(i, next_idx34) = (1-p)*(1-q);
   
 

   
    
    % Update reward matrices
    R1(i) = 1 /(state(1)+state(2));
    R2(i) = 1/ (state(1)+state(2));
    R3(i) = 1/ (state(1)+state(2));
    
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
markers = {'b*', 'ro','g^'};
colors = {'blue', 'red','green'};

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








