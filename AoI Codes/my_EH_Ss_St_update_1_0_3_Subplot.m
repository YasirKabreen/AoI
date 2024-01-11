clc;
clear;

% Ask for the values of e and q
p = input('Enter the probability p: '); % Probability of Source State p
q = input('Enter the probability q: '); % successful transmission probability q
discount = 0.95;
% Define the state space components
AoI1_values = (1:30);
B_level1 = (0:15);

% Calculate the number of states
num_states = numel(AoI1_values) * numel(B_level1);

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 2); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for b1 = B_level1
        state_space(idx, :) = [aoi1,b1];
        idx = idx + 1;
    end
end

% Create a figure
figure;

% Define markers and colors for the policy
markers = {'g^', 'ro'};
colors = {'green', 'red'};

% Iterate over different values of Probability of Source State p
e_values = [0 0.02 0.08 0.2 0.5 1];
for e_idx = 1:numel(e_values)
    e = e_values(e_idx);
    % Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
 
for i = 1:num_states
    state = state_space(i, :);
    
    % Action 1: Sending Fresh Data or Not Sending Data with Empty Update
    
    if state(2) > 0 % If battery is greater than 0
        next_state11 = [1, state(2)-1+1];
        next_state12 = [1, state(2)-1];
        next_state13 = [min(state(1) + 1, max(AoI1_values)), state(2)-1+1];
        next_state14 = [min(state(1) + 1, max(AoI1_values)), state(2)-1];
        
        % Update transition probability matrices
        P1(i, ismember(state_space, next_state11, 'rows')) = q * p *  e   ;              
        P1(i, ismember(state_space, next_state12, 'rows')) = q * p * (1-e)  ;          
        P1(i, ismember(state_space, next_state13, 'rows')) = (( 1- q * p ) *  e)   ;              
        P1(i, ismember(state_space, next_state14, 'rows')) = (( 1- q * p ) * (1-e)) ;        
        
    else % If battery is 0
        % Scenario 1: e (with energy harvesting)
        next_state11 = [min(state(1) + 1, max(AoI1_values)), min(state(2)+1, max(B_level1))];
        
        % Scenario 2: (1-e) (without energy harvesting)
        next_state12 = [min(state(1) + 1, max(AoI1_values)), state(2)];
        
        % Update transition probability matrices for Action 1 with battery = 0
        P1(i, ismember(state_space, next_state11, 'rows')) = e ;
        P1(i, ismember(state_space, next_state12, 'rows')) = 1 - e ;
    end
    
    % Action 2: Not Sending and Conserving Energy
    next_state21 = [min(state(1) + 1, max(AoI1_values)),  min(state(2)+1, max(B_level1)) ];
    next_state22 = [min(state(1) + 1, max(AoI1_values)), state(2)];
    
    % Update transition probability matrices for Action 2
    P2(i, ismember(state_space, next_state21, 'rows')) = e+P2(i, ismember(state_space, next_state21, 'rows')) ;
    P2(i, ismember(state_space, next_state22, 'rows')) = 1 - e +P2(i, ismember(state_space, next_state22, 'rows')) ;
    
    % Compute Costs for each state-action pair
    AoI_cost = state(1);
    Energy_cost_A1 = max(B_level1) - state(2);
    Energy_saved_A2 = state(2);
    w1 = 0.9; w2 = 0.1;
    w3 = 0.9; w4 = 0.1;
   if state(2) > 0 % If battery is greater than 0
    Cost_A1 = w1 * (AoI_cost -p*q*state(1))+ w2 * Energy_cost_A1;
    Cost_A2 = w3 * AoI_cost - w4 * Energy_saved_A2;
   else
    Cost_A1 = w1 * (AoI_cost )+ w2 * Energy_cost_A1;
    Cost_A2 = w3 * AoI_cost - w4 * Energy_saved_A2;
   end
    
    % Update reward matrices
    R1(i) = -Cost_A1;
    R2(i) = -Cost_A2;
   
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;


% Combine reward matrices for all sensors
R = [R1, R2];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

    
    % Create subplot
    subplot(2, 3, e_idx);
    
    % Iterate through the policy and plot the corresponding marker and color
    for i = 1:numel(policy)
        action = policy(i);
        marker = markers{action};
        color = colors{action};
        scatter(state_space(i, 2), state_space(i, 1), 50, marker, 'filled', 'MarkerEdgeColor', color);
        hold on; % To ensure points are plotted on the same graph
    end
    
    % Add labels and title
    xlabel('Battery level');
    ylabel('AoI');
    title(['EH Probability ?  = ' num2str(e)]);
    
    % Set the X and Y axis limits with a step size of 1
    axis([0 16 0 31]);
    xticks(1:1:16);
    yticks(0:1:31);
    
    % Create a legend
    if e_idx == 1
        legend('Action 1', 'Action 0');
    end
    
    % Display the plot
    grid on;
end
