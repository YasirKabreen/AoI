clc;
clear;

% Define the state space components
AoI1_values = (1:6);
AoI2_values = (1:6);
Source1_state=[1 2];
Source2_state=[1 2];



% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(Source1_state)* numel(Source2_state);

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 4); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for Ss1 = Source1_state
            for Ss2 = Source2_state
        
                        state_space(idx, :) = [aoi1, aoi2,Ss1,Ss2];
                        idx = idx + 1;
            end
        end
    end
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);
%P3 = zeros(num_states, num_states);


% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
%R3 = zeros(num_states, 1);


discount = 0.95;

% Loop through all states to calculate transitions and rewards
nn = zeros(1,num_states);
p = input('enter the probability p'); % Probability of sending
q = input('enter the probability q'); % Probability of source state
for i = 1:num_states
    state = state_space(i, :);
    
    % Sensor one transitions and rewards
  
    next_state11 = state;
    next_state12 = state;
    if(next_state11(3)==1)
    %send with probabilty p.
        next_state11(1) = 1;
        next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
        next_state11(3) = 1;
        next_state11(4) = 1;
            
         
    end
     %not send with probabilty 1-p.
     next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
     next_state12(2) = min(next_state12(2) + 1, max(AoI2_values));
    
   
    %nn = [i; next_state11];
   
    % Sensor two transitions and rewards
  next_state21 = state;
    next_state22 = state;
    
    %send with probabilty p.
    if(next_state21(4)==1)
        next_state21(2) = 1;
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
        if rand() <= 1-q % change state with probability 1-q
            next_state21(4) = 2;
            
        end
    end
     %not send with probabilty 1-p.
     next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
     next_state22(2) = min(next_state22(2) + 1, max(AoI2_values));
   
     
    
    
   
    % Find indices of next states in state_space
    next_idx11 = find(ismember(state_space, next_state11, 'rows'));
     next_idx12 = find(ismember(state_space, next_state12, 'rows'));
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));
    next_idx22 = find(ismember(state_space, next_state22, 'rows'));
    
    
    
    % Update transition probability matrices
    P1(i, next_idx11) = p;
    P1(i, next_idx12) = 1-p;
    P2(i, next_idx21) = 1-p;
    P2(i, next_idx22) = p;
   
 

   
    
    % Update reward matrices
    R1(i) = 1 / (state(1)+state(2));
    R2(i) = 1 / (state(1)+state(2));
    
    
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



% Initialize variables for simulation
current_state = [1, 1,1,1]; % Initial state (you can choose any initial state)
num_steps = 100; % Number of simulation steps


% Initialize a variable to store AoI values over time
aoi_history = zeros(num_steps, 4);

% Simulate the system
for step = 1:num_steps
    % Determine the action to take based on the current state and policy
    current_idx = find(ismember(state_space, current_state, 'rows'));
    action = policy(current_idx);
    
    % Update AoI values based on the chosen action
    if action == 1 % Sensor 1 sends
        if rand() <= p % Send with probability p
            if(current_state(3)==1)
            current_state(1) = 1;
            current_state(2) = min(current_state(2) + 1,max(AoI2_values));
            if rand() <= q % change state with probability p
            current_state(3) = 2;
             end
        else
            current_state(1) = min(current_state(1) + 1,max(AoI1_values));
             current_state(2) = min(current_state(2) + 1,max(AoI2_values));
            end
        else
             current_state(1) = min(current_state(1) + 1,max(AoI1_values));
             current_state(2) = min(current_state(2) + 1,max(AoI2_values));
        end
    elseif action == 2 % Sensor 2 sends
        if rand() <= 1-p % Send with probability p
             if(current_state(3)==1)
            current_state(2) = 1;
             current_state(1) = min(current_state(1) + 1,max(AoI1_values));
              if rand() <= 1-q % change state with probability p
            current_state(4) = 2;
             end
        else
          current_state(1) = min(current_state(1) + 1,max(AoI1_values));
           current_state(2) = min(current_state(2) + 1,max(AoI2_values));
             end
       else
          current_state(1) = min(current_state(1) + 1,max(AoI1_values));
           current_state(2) = min(current_state(2) + 1,max(AoI2_values));
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


