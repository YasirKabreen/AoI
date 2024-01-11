

starting_state_index=13;
action = policy(starting_state_index)
    reward = R(starting_state_index, action)
   
    
    % Get the probabilities for the next state given the current state and action
    prob_next_state = P(starting_state_index, :, action);
    
    % Convert the probabilities to a cumulative distribution function (CDF)
    cdf_next_state = cumsum(prob_next_state);
    
    % Generate a random number to determine the next state
    random_num = rand();
    next_state_index = find(cdf_next_state >= random_num, 1, 'first')
    
    % Update the current state for the next iteration
    current_state = StateSpace(next_state_index, :)
    
    
    
    % Plotting the Policy
figure;
plot(policy);
xlabel('State Index');
ylabel('Action');
title('Optimal Policy for AoI Minimization');
axis([0 length(StateSpace)+1 0 4]);
xticks(1:length(StateSpace));
yticks(0:4);
yticklabels({ 'No Action', 'Sensor 1', 'Sensor 2', 'Sensor 3'});
grid on;


% Create a bar plot of the policy
figure;
bar(1:num_actions, action_counts, 0.5, 'FaceColor', [0.3, 0.6, 0.9]);
xlabel('Action');
ylabel('Count');
title('Distribution of Actions in Optimal Policy');
xlim([0.5, num_actions + 0.5]);
xticks(1:num_actions);
yticks(0:max(action_counts) + 1);

% Add labels to the bars with larger numbers (in black)
for i = 1:num_actions
    text(i, action_counts(i) + 0.2, sprintf('%d', action_counts(i)), 'HorizontalAlignment', 'center', 'FontSize', 10);
end

grid on;
set(gca, 'XTickLabel', []);
box off;





figure;
subplot(2, 1, 1);
plot(policy);
xlabel('State Index');
ylabel('Action');
title('Optimal Policy for AoI Minimization');
axis([0 length(StateSpace)+1 0 4]);
xticks(1:length(StateSpace));
yticks(0:4);
yticklabels({'No Action', 'Sensor 1', 'Sensor 2', 'Sensor 3'});
grid on;

% Plotting AoI values for each sensor
subplot(2, 1, 2);
plot(1:simulation_time_steps, aoi_sensor_1, 'b', 1:simulation_time_steps, aoi_sensor_2, 'g', 1:simulation_time_steps, aoi_sensor_3, 'r');
xlabel('Time Step');
ylabel('AoI');
title('AoI Evolution for Each Sensor');
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
axis([0 simulation_time_steps 0 4]);
grid on;




target_state = [1 1 5 5 1 1];
optimal_starting_state_index = find(all(StateSpace == target_state, 2))









