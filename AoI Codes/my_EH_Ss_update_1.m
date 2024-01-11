

% The problem revolves around handling requests, considering the trade-off 
% between sending fresh data and conserving energy. We need to analyze
% whether to send data to provide users with fresh information or not to
% send, thus conserving energy and serving the user with data from a cache.

clc;
clear;
% Define the state space components
AoI1_values = (1:30);

B_level1=(0:15);



% Calculate the number of states
num_states = numel(AoI1_values) *   numel(B_level1)   ;

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


% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);



% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);



discount = 0.95;

% Loop through all states to calculate transitions and rewards
e= input('enter the probability e: '); % Probability of EH
q= input('enter the probability q: '); % Probability of Source State
for i = 1:num_states
    state = state_space(i, :);
    
  % Sensor one transitions and rewards if it Sending Fresh Data:
    next_state11 = state;
    next_state12 = state;
    next_state13 = state;
    next_state14 = state;
  
  if state(2)>=1 
        next_state11(1) = 1;
        next_state11(2)=next_state11(2);
   
       
  else
         next_state11(1) = min(next_state11(1) + 1, max(AoI1_values));
         next_state11(2)=next_state11(2)+1;
  end
  
  
    if state(2)>=1 
        next_state12(1) = 1;
        next_state12(2)=next_state12(2)-1;
   
       
        else
         next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
         
    end
    
  
      next_state13(1) = min(next_state13(1) + 1, max(AoI1_values));
      next_state13(2)=next_state13(2)+1;
    
      
      next_state14(1) = min(next_state14(1) + 1, max(AoI1_values));
        
  
 % Sensor one transitions and rewards if it is not Not Sending and Conserving Energy:
    next_state21 = state;
    next_state22 = state;
  
      
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
         next_state21(2)=min(next_state21(2)+1,max(B_level1));
    
        next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
   
    % Find indices of next states in state_space
  next_idx11 = find(ismember(state_space, next_state11, 'rows'));
       next_idx12 = find(ismember(state_space, next_state12, 'rows'));
          next_idx13 = find(ismember(state_space, next_state13, 'rows'));
             next_idx14 = find(ismember(state_space, next_state14, 'rows'));
          
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));
       next_idx22 = find(ismember(state_space, next_state22, 'rows'));
        
   
    
    % Update transition probability matrices
    P1(i, next_idx11) = e*q;
    P1(i, next_idx12) = (1-e)*q;
    P1(i, next_idx13) = e*(1-q);
    P1(i, next_idx14) = (1-e)*(1-q);
    
    P2(i, next_idx21) = e;
    P2(i, next_idx22) = 1-e;
     
   
    % Update reward matrices
    
     % Calculate AoI cost (inversely proportional to AoI)
    AoI_cost = 1 / state(1) ;

    % Define weights for AoI and Energy objectives
    w1 = 0.6; % Weight for AoI objective
    w2 = 0.4; % Weight for Energy (Battery) conservation objective

    % Calculate the energy (Battery) conservation cost
    % You can use the remaining battery level as a measure of energy conservation.
    BatteryLevel = state(2); %   state(2) represents remaining battery level
    Energy_cost = 1 - BatteryLevel / max(B_level1); % Normalize to [0, 1]

    % Calculate the overall cost using the specified weights
    Cost = w1 * AoI_cost + w2 * Energy_cost;
    R1(i) =  Cost;
    R2(i) = Cost;
    
    
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;


% Combine reward matrices for all sensors
R = [R1, R2];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);



% Create a figure
figure;

% Define markers and colors for the policy
markers = {'g^', 'ro','b*'};
colors = {'green', 'red','blue'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(policy)
    action = policy(i);
    marker = markers{action};
    color = colors{action};
    scatter(state_space(i, 2), state_space(i, 1), 50, marker, 'filled', 'MarkerEdgeColor', color);
    hold on; % To ensure points are plotted on the same graph
end

% Add labels
xlabel('Battery level');
ylabel('AoI');
% Set the X and Y axis limits with a step size of 1
axis([1 16 0 31]);
xticks(1:1:35);
yticks(0:1:35);
% Create a legend
legend('Action 1', 'Action 2');

% Display the plot
grid on;








