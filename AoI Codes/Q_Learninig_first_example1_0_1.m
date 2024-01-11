 

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
         next_state11(2)=next_state11(2);
  end
  
  
 if state(2)>=1 
        next_state12(1) = 1;
        next_state12(2)=max(next_state12(2)-1,0);
   
       
        else
         next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
          next_state12(2)=max(next_state12(2)-1,0);
         
 end
 
  
              next_state13(1) = min(next_state13(1) + 1, max(AoI1_values));
   
       
              next_state14(1) = min(next_state14(1) + 1, max(AoI1_values));
              next_state14(2)=max(next_state14(2)-1,0);
      
        
  
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
    P1(i, next_idx11) = P1(i, next_idx11)+ e*q;
    P1(i, next_idx12) = P1(i, next_idx12)+(1-e)*q;
    P1(i, next_idx13) = P1(i, next_idx13)+e*(1-q);
    P1(i, next_idx14) = P1(i, next_idx14)+(1-e)*(1-q);
    
    P2(i, next_idx21) = P2(i, next_idx21)+e;
    P2(i, next_idx22) = P2(i, next_idx22)+1-e;
     
   
    % Update reward matrices
    
     % Calculate AoI cost (inversely proportional to AoI)
    AoI_cost = 1 / state(1) ;

    % Define weights for AoI and Energy objectives
    w1 = 0.1; % Weight for AoI objective
    w2 = 0.9; % Weight for Energy (Battery) conservation objective

    % Calculate the energy (Battery) conservation cost
    % You can use the remaining battery level as a measure of energy conservation.
    BatteryLevel = state(2); %   state(2) represents remaining battery level
    Energy_cost =  max(B_level1) - BatteryLevel;
   

    % Calculate the overall cost using the specified weights
    Cost =  .9*state(1)+.1*Energy_cost;
    R1(i) = 1/Cost;
    R2(i) = 1/Cost;
    
     
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;


% Combine reward matrices for all sensors
R = [R1, R2];




% Q-Learning parameters
alpha = 0.1; % learning rate
gamma = 0.9 ; % discount factor
epsilon = 0.4; % exploration rate
num_episodes = 1000; % number of episodes for training

% Initialize Q-table with zeros
Q = zeros(num_states, 2); % Two actions: Send Fresh Data or Conserve Energy

for episode = 1:num_episodes
    % Assuming starting state is random
    state_idx = randi([1, num_states]);
    
    done = false; % to denote if episode has ended
    while ~done
        % Epsilon-greedy action selection
        if rand() < epsilon
            action = randi([1, 2]); % choose random action
        else
            [~, action] = max(Q(state_idx, :)); % choose best action based on Q-values
        end
        
        % Given the selected action, get the next state using the transition probabilities
        next_state_probs = P(state_idx, :, action);
        next_state_idx = find(mnrnd(1, next_state_probs) == 1); % choose a next state
        
        % Get the immediate reward for the selected action
        reward = R(state_idx, action);
        
        % Q-learning update rule
        best_future_value = max(Q(next_state_idx, :));
        Q(state_idx, action) = Q(state_idx, action) + alpha * (reward + gamma * best_future_value - Q(state_idx, action));
        
        % Move to the next state
        state_idx = next_state_idx;
        
        % Termination condition (you can use a fixed length for episodes or define another criterion)
        if rand() > 0.95 % 5% chance to end the episode, adapt this as needed
            done = true;
        end
    end
end

% Derive the optimal policy from the Q-table
[~, optimal_policy] = max(Q, [], 2);

% Create a figure
figure;

% Define markers and colors for the policy
markers = {'g^', 'ro','b*'};
colors = {'green', 'red','blue'};

% Iterate through the policy and plot the corresponding marker and color
for i = 1:numel(optimal_policy)
    action = optimal_policy(i);
    marker = markers{action};
    color = colors{action};
    scatter(state_space(i, 2), state_space(i, 1), 50, marker, 'filled', 'MarkerEdgeColor', color);
    hold on; % To ensure points are plotted on the same graph
end

% Add labels
xlabel('Battery level');
ylabel('AoI');
% Set the X and Y axis limits with a step size of 1
axis([0 16 0 31]);
xticks(1:1:35);
yticks(0:1:35);
% Create a legend
legend('Action 1', 'Action 2');

% Display the plot
grid on;








