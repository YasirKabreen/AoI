clc;
clear;
% Define the state space components
AoI1_values = (1:10);

B_level1=(0:6);
B_level2=(0:6);


% Calculate the number of states
num_states = numel(AoI1_values) *   numel(B_level1) * numel(B_level2) ;

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 3); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
     
        for b1 = B_level1
            for b2 = B_level2
                        state_space(idx, :) = [aoi1,b1,b2];
                        idx = idx + 1;
            end
        end
     
end
state_space

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

q = 0.6;
e=.7;
for i = 1:num_states
    state = state_space(i, :);
    
% when send Sensor one transitions and rewards
    next_state11 = state;
    next_state12 = state;
    next_state13 = state;
    next_state14 = state;
  %with probability e*q
  if state(2)>=1 
        next_state11(1) = 1;
        next_state11(2)=next_state11(2);
        next_state11(3)=min(next_state11(3)+1,max(B_level2));
   
       
  else
         next_state11(1) = min(next_state11(1) + 1, max(AoI1_values));
         next_state11(2)=min(next_state11(2)+1,max(B_level1));
         next_state11(3)=min(next_state11(3)+1,max(B_level2));
  end
  
    %with probability 1-e*q
    if state(2)>=1 
        next_state12(1) = 1;
        next_state12(2) =  next_state12(2)-1;
        next_state12(3) =  next_state11(3);
        
   
       
        else
         next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
         
    end
    
      %with probability e*1-q
      next_state13(1) = min(next_state13(1) + 1, max(AoI1_values));
      next_state13(2)=min(next_state13(2)+1, max(B_level1));
      next_state13(3)=min(next_state13(3)+1, max(B_level2));
    
      %with probability 1-e*1-q
      next_state14(1) = min(next_state14(1) + 1, max(AoI1_values));
        
  
 %when sending  Sensor two transitions and rewards
        next_state21 = state;
        next_state22 = state;
        next_state23 = state;
        next_state24 = state;
  %with probability e*1-q
  if state(3)>=1 
        next_state21(1) = 1;
        next_state21(3) = next_state21(3);
        next_state21(2) = min(next_state21(2)+1,max(B_level1));
   
       
  else
         next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
         next_state21(2) = min(next_state21(2)+1,max(B_level1));
         next_state21(3) = min(next_state21(3)+1,max(B_level2));
  end
  
    %with probability 1-e*1-q
    if state(3)>=1 
        next_state22(1) = 1;
        next_state22(3) =  next_state22(3)-1;
        next_state22(2) =  next_state21(2);
        
   
       
        else
         next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
         
    end
    
      %with probability e*q
      next_state23(1) = min(next_state23(1) + 1, max(AoI1_values));
      next_state23(2) = min(next_state23(2)+1, max(B_level1));
      next_state23(3) = min(next_state23(3)+1, max(B_level2));
    
      %with probability 1-e*q
      next_state24(1) = min(next_state24(1) + 1, max(AoI1_values));
   
 %when not requested to send from the BS
        next_state31 = state;
        next_state32 = state;
     %with probability e
      next_state31(1) = min(next_state31(1) + 1, max(AoI1_values));
      next_state31(2) = min(next_state31(2)+1, max(B_level1));
      next_state31(3) = min(next_state31(3)+1, max(B_level2));
      
      %with probability 1-e
      next_state32(1) = min(next_state32(1) + 1, max(AoI1_values));
      
      
      
      
      
      
    % % Find indices of next states in state_space
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
    
    % Update transition probability matrices
    P1(i, next_idx11) = e*q;
    P1(i, next_idx12) = (1-e)*q;
    P1(i, next_idx13) = e*(1-q);
    P1(i, next_idx14) = (1-e)*(1-q);
    
    P2(i, next_idx21) =  e*(1-q);
    P2(i, next_idx22) = (1-e)*(1-q);
    P2(i, next_idx23) = e*q;
    P2(i, next_idx24) = (1-e)*q;
        
    P3(i, next_idx31) = e ;
    P3(i, next_idx32) = (1-e) ;
   
    
      % Update reward matrices
    bb12=(next_state11(2)+next_state12(2)+next_state13(2)+next_state14(2))/4;
     bb13=(next_state11(3)+next_state12(3)+next_state13(3)+next_state14(3))/4;
     
      bb22=(next_state21(2)+next_state22(2)+next_state23(2)+next_state24(2))/4;
     bb23=(next_state21(3)+next_state22(3)+next_state23(3)+next_state24(3))/4;
     
     bb32=(next_state31(2)+next_state32(2))/2;
     bb33=(next_state31(3)+next_state32(3))/2;
     
    bb2=(next_state21(2)+next_state22(2))/2;
    R1(i) =  (bb12+bb13) /(next_state11(1)+next_state12(1)+next_state13(1)+next_state14(1))/4;
    R2(i) = (bb22+bb23)/ (next_state21(1)+next_state22(1)+next_state23(1)+next_state24(1))/4;
    R3(i) = (bb32+bb33)/ (next_state31(1)+next_state32(1))/2;
    
    
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
    scatter(i, state_space(i, 1), 100, marker, 'filled', 'MarkerEdgeColor', color);
    hold on; % To ensure points are plotted on the same graph
end

% Add labels
xlabel('AoI1');
ylabel('AoI2');
% Set the X and Y axis limits with a step size of 1
axis([0 210 0 21]);
xticks(0:1:21);
yticks(0:1:21);
% Create a legend
legend('Action 1', 'Action 2');

% Display the plot
grid on;








