clc;
clear;

% Define the state space components
AoI1_values = 1:5;
AoI2_values = 1:5;
Battery_levels = 0:5;  % Battery levels including 0
e = 0.2;  % Example energy harvesting probability
discount=0.9;
% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(Battery_levels)^3;

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 5);  % Including battery levels for three sensors

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for bat1 = Battery_levels
            for bat2 = Battery_levels
                for bat3 = Battery_levels
                    state_space(idx, :) = [aoi1, aoi2, bat1, bat2, bat3];
                    idx = idx + 1;
                end
            end
        end
    end
end

% Initialize transition probability matrices
num_actions = 3;
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);
P3 = zeros(num_states, num_states);
P4 = zeros(num_states, num_states);

% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
R3 = zeros(num_states, 1);
R4 = zeros(num_states, 1);
q = input('Enter the probability q: ');

for i = 1:num_states
   state = state_space(i, :);
    
    % Action 1: Sending Fresh Data from sensor 1
    
    if state(3) > 0 % If battery is greater than 0
        next_state11 = [1,min(state(2) + 1, max(AoI2_values)), state(3)-1+1,min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
        next_state12 = [1,min(state(2) + 1, max(AoI2_values)), state(3)-1,min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
        next_state13 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
        next_state14 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
        
        % Update transition probability matrices
        P1(i, ismember(state_space, next_state11, 'rows')) = e * q +  P1(i, ismember(state_space, next_state11, 'rows'));
        P1(i, ismember(state_space, next_state12, 'rows')) = (1-e) * q +P1(i, ismember(state_space, next_state12, 'rows'));
        P1(i, ismember(state_space, next_state13, 'rows')) = e * (1-q) +P1(i, ismember(state_space, next_state13, 'rows'));
        P1(i, ismember(state_space, next_state14, 'rows')) = (1-e) * (1-q) + P1(i, ismember(state_space, next_state14, 'rows'));
        
    else % If battery is 0
        % Scenario 1: e (with energy harvesting)
         next_state11 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
       
        
        % Scenario 2: (1-e) (without energy harvesting)
         next_state12 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
        
        % Update transition probability matrices for Action 1 with battery = 0
        P1(i, ismember(state_space, next_state11, 'rows')) = e+ P1(i, ismember(state_space, next_state11, 'rows'));
        P1(i, ismember(state_space, next_state12, 'rows')) = 1 - e+ P1(i, ismember(state_space, next_state12, 'rows'));
    end
    
    % Action 2: Sending Fresh Data from sensor 2
    if state(4) > 0 % If battery is greater than 0
   next_state21 = [min(state(1) + 1, max(AoI1_values)),1, min(state(3)+1, max(Battery_levels)),state(4)-1+1,min(state(5)+1, max(Battery_levels))];
    next_state22 = [min(state(1) + 1, max(AoI1_values)),1, min(state(3), max(Battery_levels)),state(4)-1,min(state(5), max(Battery_levels))];
      next_state23 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
        next_state24 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
      
          % Update transition probability matrices
        P2(i, ismember(state_space, next_state21, 'rows')) = e * q +  P2(i, ismember(state_space, next_state21, 'rows'));
        P2(i, ismember(state_space, next_state22, 'rows')) = (1-e) * q +P2(i, ismember(state_space, next_state22, 'rows'));
        P2(i, ismember(state_space, next_state23, 'rows')) = e * (1-q) +P2(i, ismember(state_space, next_state23, 'rows'));
        P2(i, ismember(state_space, next_state24, 'rows')) = (1-e) * (1-q) + P2(i, ismember(state_space, next_state24, 'rows'));
   else % If battery is 0
        % Scenario 1: e (with energy harvesting)
         next_state21 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
       
        
        % Scenario 2: (1-e) (without energy harvesting)
         next_state22 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
        
        % Update transition probability matrices for Action 1 with battery = 0
        P2(i, ismember(state_space, next_state21, 'rows')) = e+ P2(i, ismember(state_space, next_state21, 'rows'));
        P2(i, ismember(state_space, next_state22, 'rows')) = 1 - e+ P2(i, ismember(state_space, next_state22, 'rows'));
    end
    
    
 
%Action 3: Sending Fresh Data from both sensor 1 and 2

    if state(5) > 0 % If battery is greater than 0
   next_state31 = [1,1, min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),state(5)-1+1];
    next_state32 = [1,1, min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),state(5)-1];
      next_state33 = [1,min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),state(5)-1+1];
        next_state34 = [1,min(state(2) + 1, max(AoI2_values)), min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),state(5)-1];
       next_state35 = [min(state(1) + 1, max(AoI1_values)),1, min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),state(5)-1+1];
    next_state36 = [min(state(1) + 1, max(AoI1_values)),1, min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),state(5)-1];
      next_state37 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
        next_state38 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(4), max(Battery_levels)),min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
        
          % Update transition probability matrices
        P3(i, ismember(state_space, next_state31, 'rows')) = (1-q) * (1-q)*e       + P3(i, ismember(state_space, next_state31, 'rows'));
        P3(i, ismember(state_space, next_state32, 'rows')) = (1-q) * (1-q)*(1-e)   + P3(i, ismember(state_space, next_state32, 'rows'));
        P3(i, ismember(state_space, next_state33, 'rows')) = (1-q) * q*e           + P3(i, ismember(state_space, next_state33, 'rows'));
        P3(i, ismember(state_space, next_state34, 'rows')) = (1-q) * q*(1-e)       + P3(i, ismember(state_space, next_state34, 'rows'));
      
        P3(i, ismember(state_space, next_state35, 'rows')) = q * (1-q)*e                 +  P3(i, ismember(state_space, next_state35, 'rows'));
        P3(i, ismember(state_space, next_state36, 'rows')) = q * (1-q)*(1-e)             +  P3(i, ismember(state_space, next_state36, 'rows'));
        P3(i, ismember(state_space, next_state37, 'rows')) = q * q*e              +  P3(i, ismember(state_space, next_state37, 'rows'));
        P3(i, ismember(state_space, next_state38, 'rows')) = q * q  *(1-e)       +  P3(i, ismember(state_space, next_state38, 'rows'));
   else % If battery is 0
        % Scenario 1: e (with energy harvesting)
         next_state31 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
       
        
        % Scenario 2: (1-e) (without energy harvesting)
         next_state32 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3), max(Battery_levels)),min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
        
        % Update transition probability matrices for Action 1 with battery = 0
        P3(i, ismember(state_space, next_state31, 'rows')) = e+ P3(i, ismember(state_space, next_state31, 'rows'));
        P3(i, ismember(state_space, next_state32, 'rows')) = 1 - e+ P3(i, ismember(state_space, next_state32, 'rows'));
    end
    
    
%Action 4: Not sending any data (conserving energy)
     % Scenario 1: e (with energy harvesting)
         next_state41 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(3)+1, max(Battery_levels)),min(state(4)+1, max(Battery_levels)),min(state(5)+1, max(Battery_levels))];
       
        
        % Scenario 2: (1-e) (without energy harvesting)
         next_state42 = [min(state(1) + 1, max(AoI1_values)),min(state(2) + 1, max(AoI2_values)), min(state(4), max(Battery_levels)),min(state(4), max(Battery_levels)),min(state(5), max(Battery_levels))];
        
        % Update transition probability matrices for Action 1 with battery = 0
        P4(i, ismember(state_space, next_state41, 'rows')) = e+ P4(i, ismember(state_space, next_state41, 'rows'));
        P4(i, ismember(state_space, next_state42, 'rows')) = 1 - e+ P4(i, ismember(state_space, next_state42, 'rows'));
    
    % Compute Costs for each state-action pair
    AoI_cost1 = state(1);
    AoI_cost2 = state(2);
    Energy_cost_A1 = max(Battery_levels) - state(3);% Define the cost for using energy
    Energy_cost_A2 = max(Battery_levels) - state(4);% Define the cost for using energy
    Energy_cost_A3 = max(Battery_levels) - state(5);% Define the cost for using energy
    Energy_saved_A4 = state(3)+state(4)+state(5);   % Define the cost for saving energy
    w1 = 0.9; w2 = 0.1;
    w3 = 0.5; w4 = 0.1;
    Cost_A1 = w2 * AoI_cost1 + w2 * AoI_cost2 + w4 * Energy_cost_A1;
    Cost_A2 = w1 * AoI_cost1 + w2 * AoI_cost2 + w4 * Energy_cost_A2;
    Cost_A3 = w3 * AoI_cost1 + w3 * AoI_cost2 + w4 * Energy_cost_A3;
    Cost_A4 = w3 * AoI_cost1 + w3 * AoI_cost2 - w4 * Energy_saved_A4;
    
    % Update reward matrices
    R1(i) = -Cost_A1;
    R2(i) = -Cost_A2;
    R3(i) = -Cost_A3;
    R4(i) = -Cost_A4;
        
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;
P(:,:,4) = P4;


% Combine reward matrices for all sensors
R = [R1, R2,R3,R4];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);


% ...

% Create a figure
figure;

% Define markers and colors for the policy
markers = {'ro','b*', 'g^'};
colors = {'red','blue', 'green'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(policy)
    action = policy(i);
    marker = markers{action};
    color = colors{action};

    % Check battery levels
    if state_space(i, 3) == 1 && state_space(i, 4) == 5 && state_space(i, 5) == 1
        scatter(state_space(i, 1), state_space(i, 2), 50, marker, 'filled', 'MarkerEdgeColor', color);
        hold on; % To ensure points are plotted on the same graph
    end
end

% Add labels
xlabel('AoI1');
ylabel('AoI2');
% Set the X and Y axis limits with a step size of 1
axis([0 6 0 6]);
xticks(0:1:6);
yticks(0:1:6);
% Create a legend
legend('Action 1', 'Action 2','Action 3');

% Display the plot
grid on;

% Create a figure
figure;

% Define markers and colors for the policy
markers = {'ro','b*', 'g^'};
colors = {'red','blue', 'green'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(policy)
    action = policy(i);
    marker = markers{action};
    color = colors{action};
    scatter(state_space(i, 1), state_space(i, 2), 50, marker, 'filled', 'MarkerEdgeColor', color);
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
legend('Action 1', 'Action 2','Action 3');

% Display the plot
grid on;



