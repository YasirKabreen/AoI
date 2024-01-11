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


% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);



discount = 0.95;

% Loop through all states to calculate transitions and rewards
nn = zeros(1,num_states);
p=0.9;
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
    
   
    %nn = [i; next_state11];
   
    % Sensor two transitions and rewards
  next_state21 = state;
    next_state22 = state;
    
    %send with probabilty p.
        next_state21(2) = 1;
        next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
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
    P2(i, next_idx21) = p;
    P2(i, next_idx22) = 1-p;
 

   
    
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

% Initialize the total AoI
total_AoI = 0;

% Specify the number of time steps for the simulation
num_time_steps = 100; 

% Specify the initial state 
initial_state = [1 1]; 

% Simulate the system
current_state = initial_state;
vv = [];
for t = 1:num_time_steps
    % Use the policy to determine the action to take
    action = policy(find(all(bsxfun(@eq, state_space, current_state), 2), 1));

    
    % Transition to the next state based on the chosen action
    if (action == 1)
        % Perform Sensor 1 action
        if (current_state(4) == 1 || current_state(4) == 2)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        if (current_state(5) == 1 || current_state(5) == 2)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI1_values));
        end
        if (current_state(6) == 1 || current_state(6) == 2)
            current_state(3) = 1;
        else
            current_state(3) = min(current_state(3) + 1, max(AoI1_values));
        end
        
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
        current_state(6) = mod(current_state(6) + 1, num_zones + 1);
        current_state(6) = max(current_state(6), 1); % Ensure it's at least 1
    elseif (action == 2)
        % Perform Sensor 2 action
        if (current_state(4) == 3)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        if (current_state(5) == 3)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI1_values));
        end
        if (current_state(6) == 3)
            current_state(3) = 1;
        else
            current_state(3) = min(current_state(3) + 1, max(AoI1_values));
        end
        
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
        current_state(6) = mod(current_state(6) + 1, num_zones + 1);
        current_state(6) = max(current_state(6), 1); % Ensure it's at least 1
    elseif (action == 3)
        % Perform Sensor 3 action
        if (current_state(4) == 5)
            current_state(1) = 1;
        else
            current_state(1) = min(current_state(1) + 1, max(AoI1_values));
        end
        if (current_state(5) == 5)
            current_state(2) = 1;
        else
            current_state(2) = min(current_state(2) + 1, max(AoI1_values));
        end
        if (current_state(6) == 5)
            current_state(3) = 1;
        else
            current_state(3) = min(current_state(3) + 1, max(AoI1_values));
        end
        
        current_state(4) = mod(current_state(4) + 1, num_zones + 1);
        current_state(4) = max(current_state(4), 1); % Ensure it's at least 1
        current_state(5) = mod(current_state(5) + 1, num_zones + 1);
        current_state(5) = max(current_state(5), 1); % Ensure it's at least 1
        current_state(6) = mod(current_state(6) + 1, num_zones + 1);
        current_state(6) = max(current_state(6), 1); % Ensure it's at least 1
    end
    
    % Calculate and update the total AoI
    total_AoI = total_AoI + current_state(1) + current_state(2) + current_state(3); % Assuming AoI1, AoI2, and AoI3 are stored in current_state
    
   
    
    vv = [vv; current_state];
end

fprintf('Total AoI: %d\n', total_AoI);
