clc;
clear;

% Define the state space components
AoI1_values = (1:30);
B_level1=(0:15);

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
not_sent=[];
sent=[];
% Loop through the EH probabilities from 0.05 to 1 with a step of 0.2
for e = 0.005:0.005:1

    % Reset transition probability matrices for each loop
    P1 = zeros(num_states, num_states);
    P2 = zeros(num_states, num_states);

    % Reset reward matrix for each loop
    R1 = zeros(num_states, 1);
    R2 = zeros(num_states, 1);

    discount = 0.95;
    q = 0.9; % Assuming a constant value for q based on your provided code

   for i = 1:num_states
    state = state_space(i, :);
    
    % Action 1: Sending Fresh Data or Not Sending Data with Empty Update
    
    if state(2) > 0 % If battery is greater than 0
        next_state11 = [1, state(2)-1+1];
        next_state12 = [1, state(2)-1];
        next_state13 = [min(state(1) + 1, max(AoI1_values)), state(2)-1+1];
        next_state14 = [min(state(1) + 1, max(AoI1_values)), state(2)-1];
        
        % Update transition probability matrices
        P1(i, ismember(state_space, next_state11, 'rows')) = e * q+  P1(i, ismember(state_space, next_state11, 'rows'));
        P1(i, ismember(state_space, next_state12, 'rows')) = (1-e) * q+P1(i, ismember(state_space, next_state12, 'rows'));
        P1(i, ismember(state_space, next_state13, 'rows')) = e * (1-q)+P1(i, ismember(state_space, next_state13, 'rows'));
        P1(i, ismember(state_space, next_state14, 'rows')) = (1-e) * (1-q)+ P1(i, ismember(state_space, next_state14, 'rows'));
        
    else % If battery is 0
        % Scenario 1: e (with energy harvesting)
        next_state11 = [min(state(1) + 1, max(AoI1_values)), min(state(2)+1, max(B_level1))];
        
        % Scenario 2: (1-e) (without energy harvesting)
        next_state12 = [min(state(1) + 1, max(AoI1_values)), state(2)];
        
        % Update transition probability matrices for Action 1 with battery = 0
        P1(i, ismember(state_space, next_state11, 'rows')) = e+ P1(i, ismember(state_space, next_state11, 'rows'));
        P1(i, ismember(state_space, next_state12, 'rows')) = 1 - e+ P1(i, ismember(state_space, next_state12, 'rows'));
     end
   
    
    % Action 2: Not Sending and Conserving Energy
    next_state21 = [min(state(1) + 1, max(AoI1_values)),  min(state(2)+1, max(B_level1)) ];
    next_state22 = [min(state(1) + 1, max(AoI1_values)), state(2)];
    
    % Update transition probability matrices for Action 2
    P2(i, ismember(state_space, next_state21, 'rows')) = e+   P2(i, ismember(state_space, next_state21, 'rows'));
    P2(i, ismember(state_space, next_state22, 'rows')) = 1 - e+P2(i, ismember(state_space, next_state22, 'rows'));
    
    % Compute Costs for each state-action pair
    AoI_cost = state(1);
    Energy_cost_A1 = max(B_level1) - state(2);
    Energy_saved_A2 = state(2);
    w1 = 0.9; w2 = 0.1;
    w3 = 0.9; w4 = 0.1;
    Cost_A1 = w1 * AoI_cost + w2 * Energy_cost_A1;
    Cost_A2 = w3 * AoI_cost - w4 * Energy_saved_A2;
    
    % Update reward matrices
    R1(i) = -Cost_A1;
    R2(i) = -Cost_A2;
    
    % ... [your existing transition update code here]
    end

    

    % Combine transition probability matrices for all sensors
    P(:,:,1) = P1;
    P(:,:,2) = P2;

    % Combine reward matrices for all sensors
    R = [R1, R2];

    % Perform value iteration to find the optimal value function and policy
    [V, policy] = mdp_policy_iteration(P, R, discount);

    % Calculate the number of states that are not sent based on your policy
    not_sent_states = sum(policy == 2); % Assuming 2 indicates not sent
      sent_states = sum(policy == 1); % Assuming 1 indicates   sent
    fprintf('For e = %.3f, number of states not sent = %d\n', e, not_sent_states);
    not_sent=[not_sent,not_sent_states];
    sent=[sent,sent_states];
end

% Data
e_values = 0.005:0.005:1;
states_not_sent = not_sent;
states_sent = 480-not_sent;
figure
subplot(2,2,1)
plot(e_values, states_not_sent, '-');
axis([0 1 0 480]);
xlabel('Probability of EH (e)');
ylabel('Number of States Not Sent');
title('Relationship between Probability of EH and States Not Sent');
grid on;

subplot(2,2,2)
bar(e_values, states_not_sent);
axis([0 1 0 480]);
xlabel('Probability of EH (e)');
ylabel('Number of States Not Sent');
title('Relationship between Probability of EH and States Not Sent');
grid on;
subplot(2,2,3)
plot(e_values, states_sent, '-');
axis([0 1 0 480]);
xlabel('Probability of EH (e)');
ylabel('Number of States Sent');
title('Relationship between Probability of EH and States Sent');
grid on;

subplot(2,2,4)
bar(e_values, states_sent);
axis([0 1 0 480]);
xlabel('Probability of EH (e)');
ylabel('Number of States Sent');
title('Relationship between Probability of EH and States Sent');
grid on;

% Create a new figure for the modified scatter plot
figure;

% Loop through each e value
for idx = 1:length(e_values)
    e = e_values(idx);
    
    % Scatter plot for states not sent
    y_not_sent = 1:states_not_sent(idx);
    x_not_sent = e * ones(size(y_not_sent));
    scatter(x_not_sent, y_not_sent, 15,'ro', 'filled');
    hold on;

    % Scatter plot for states sent
    y_sent = (states_not_sent(idx) + 1):480;
    x_sent = e * ones(size(y_sent));
    scatter(x_sent, y_sent,15, 'g^', 'filled');
end

% Add labels, title, and legend
axis([0 1 50 110]);
xlabel('Probability of EH (e)');
ylabel('Number of States');
title('Relationship between Probability of EH and States Sent/Not Sent');
legend('States Not Sent', 'States Sent');

% Display the plot with grid
grid on;


% Data
e_values = 0.05:0.02:0.99;
states_not_sent = [90, 84, 81, 76, 74, 72, 69, 68, 66, 64, 63, 63, 62, 61, 61, 61, 60, 60, 60, 60, 60, 59, 58, 58, 58, 58, 58, 58, 58, 58, 57, 57, 57, 57, 57, 56, 56, 56, 56, 55, 55, 55, 55, 54, 54, 54, 54, 53];

%{ 
Interpolation
estimated_states_not_sent = interp1(e_values, states_not_sent, 0.7);

fprintf('Estimated number of states not sent for e = 0.9 is: %f\n', estimated_states_not_sent);
%}




