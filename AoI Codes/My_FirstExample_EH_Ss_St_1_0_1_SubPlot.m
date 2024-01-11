clc;
clear;

% Define the state space components
AoI1_values = 1:10;
AoI2_values = 1:10;
B = 3; % Maximum battery level
Battery_values = 0:B;

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(Battery_values)^3;

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 5); 

% Generate all combinations of states
idx = 1;
for b1 = Battery_values
    for b2 = Battery_values
        for b3 = Battery_values
            for aoi1 = AoI1_values
                for aoi2 = AoI2_values
                    state_space(idx, :) = [aoi1, aoi2, b1, b2, b3];
                    idx = idx + 1;
                end
            end
        end
    end
end

% Initialize transition probability matrices and reward matrix
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);
P3 = zeros(num_states, num_states);

R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
R3 = zeros(num_states, 1);
discount=0.95;
p = input('Enter the probability p: ');
q = 1;
p_eh = 0.2; % Probability of energy harvesting

empty_penalty = 1e3; % High penalty for empty battery

for i = 1:num_states
    state = state_space(i, :);
    delta1 = state(1);
    delta2 = state(2);
    b1 = state(3);
    b2 = state(4);
    b3 = state(5);

    % Battery penalty per sensor
if b1 == 0
    battery_penalty1 = empty_penalty;
else
    battery_penalty1 = b1 * 0.05;
end

if b2 == 0
    battery_penalty2 = empty_penalty;
else
    battery_penalty2 = b2 * 0.05;
end

if b3 == 0
    battery_penalty3 = empty_penalty;
else
    battery_penalty3 = b3 * 0.05;
end


    % Cost computation with battery penalty
    Cost1 = (delta1 + delta2)/2 + 1 - 0.5 * p * q * delta1 + battery_penalty1;
    Cost2 = (delta1 + delta2)/2 + 1 - 0.5 * p * q * delta2 + battery_penalty2;
    Cost3 = (delta1 + delta2)/2 + 1 - 0.5 * (1-p) * q * (delta1 + delta2) + battery_penalty3;

    % Rewards
    R1(i) = -Cost1;
    R2(i) = -Cost2;
    R3(i) = -Cost3;

   % New states after considering energy harvesting for each action:

% Define the maximum AoI values for clarity
max_AoI1 = max(AoI1_values);
max_AoI2 = max(AoI2_values);

% For successful transmissions with energy harvesting:
new_state1_eh = [1, min(delta2+1, max_AoI2), min(max(b1-1, 0)+1, B), min(b2+1, B), min(b3+1, B)];
new_state2_eh = [min(delta1+1, max_AoI1), 1, min(b1+1, B), min(max(b2-1, 0)+1, B), min(b3+1, B)];
new_state3_eh = [1, 1, min(b1+1, B), min(b2+1, B), min(max(b3-1, 0)+1, B)];

% For unsuccessful transmissions with energy harvesting:
new_state_fail1_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), min(max(b1-1, 0)+1, B), min(b2+1, B), min(b3+1, B)];
new_state_fail2_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), min(b1+1, B), min(max(b2-1, 0)+1, B), min(b3+1, B)];
new_state_fail3_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), min(b1+1, B), min(b2+1, B), min(max(b3-1, 0)+1, B)];

% For successful transmissions without energy harvesting:
new_state1_no_eh = [1, min(delta2+1, max_AoI2), max(b1-1, 0), b2, b3];
new_state2_no_eh = [min(delta1+1, max_AoI1), 1, b1, max(b2-1, 0), b3];
new_state3_no_eh = [1, 1, b1, b2, max(b3-1, 0)];

% For unsuccessful transmissions without energy harvesting:
new_state_fail1_no_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), max(b1-1, 0), b2, b3];
new_state_fail2_no_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), b1, max(b2-1, 0), b3];
new_state_fail3_no_eh = [min(delta1+1, max_AoI1), min(delta2+1, max_AoI2), b1, b2, max(b3-1, 0)];


% Finding the index for each new state in the state space after EH:
idx1_eh = find(ismember(state_space, new_state1_eh, 'rows'));
idx2_eh = find(ismember(state_space, new_state2_eh, 'rows'));
idx3_eh = find(ismember(state_space, new_state3_eh, 'rows'));

idx_fail1_eh = find(ismember(state_space, new_state_fail1_eh, 'rows'));
idx_fail2_eh = find(ismember(state_space, new_state_fail2_eh, 'rows'));
idx_fail3_eh = find(ismember(state_space, new_state_fail3_eh, 'rows'));

% Finding the index for each new state in the state space without EH:
idx1_no_eh = find(ismember(state_space, new_state1_no_eh, 'rows'));
idx2_no_eh = find(ismember(state_space, new_state2_no_eh, 'rows'));
idx3_no_eh = find(ismember(state_space, new_state3_no_eh, 'rows'));

idx_fail1_no_eh = find(ismember(state_space, new_state_fail1_no_eh, 'rows'));
idx_fail2_no_eh = find(ismember(state_space, new_state_fail2_no_eh, 'rows'));
idx_fail3_no_eh = find(ismember(state_space, new_state_fail3_no_eh, 'rows'));


% Updating the transition probabilities for each state considering EH:
P1(i, idx1_eh) = p * q * p_eh;
P1(i, idx_fail1_eh) = (1 - p*q) * p_eh;

P2(i, idx2_eh) = p * q * p_eh;
P2(i, idx_fail2_eh) = (1 - p*q) * p_eh;

P3(i, idx3_eh) = (1-p) * q * p_eh;
P3(i, idx_fail3_eh) = (1 - (1-p)*q) * p_eh;

% Updating the transition probabilities for each state considering no EH:
P1(i, idx1_no_eh) = p * q * (1 - p_eh);
P1(i, idx_fail1_no_eh) = (1 - p*q) * (1 - p_eh);

P2(i, idx2_no_eh) = p * q * (1 - p_eh);
P2(i, idx_fail2_no_eh) = (1 - p*q) * (1 - p_eh);

P3(i, idx3_no_eh) = (1-p) * q * (1 - p_eh);
P3(i, idx_fail3_no_eh) = (1 - (1-p)*q) * (1 - p_eh);


end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;

% Combine reward matrices for all sensors
R = [R1, R2, R3];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

figure;

% Define battery level combinations to cover major changes
battery_combinations = [
    1,1,1;
    0,0,2;
    0,2,2;
    2,2,0;
    2,0,2;
    2,2,2;
    0,2,3;
    0,3,2;
    2,0,3;
    3,0,2;
    2,3,0;
    3,2,0;
    1,3,3;
    3,1,3;
    3,3,1;
    3,3,3;
];

for k = 1:size(battery_combinations, 1)
    b_combination = battery_combinations(k, :);

    % Find indices corresponding to the battery combination
    selected_indices = find(...
        state_space(:,3) == b_combination(1) & ...
        state_space(:,4) == b_combination(2) & ...
        state_space(:,5) == b_combination(3));

    % Extract AoI1, AoI2, and corresponding policy choices
    AoI1_selected = state_space(selected_indices, 1);
    AoI2_selected = state_space(selected_indices, 2);
    policy_selected = policy(selected_indices);

    % Create the scatter subplot
    subplot(4,4,k);
    hold on;
    
    % Scatter for each action
    scatter(AoI1_selected(policy_selected == 1), AoI2_selected(policy_selected == 1), 'r', 'filled'); % Action 1 in red
    scatter(AoI1_selected(policy_selected == 2), AoI2_selected(policy_selected == 2), 'b', 'filled'); % Action 2 in blue
    scatter(AoI1_selected(policy_selected == 3), AoI2_selected(policy_selected == 3), 'g', 'filled'); % Action 3 in green
    
    xlabel('AoI1');
    ylabel('AoI2');
    title(sprintf('b1=%d, b2=%d, b3=%d', b_combination(1), b_combination(2), b_combination(3)));
    xlim([min(AoI1_values), max(AoI1_values)]);
    ylim([min(AoI2_values), max(AoI2_values)]);
    axis([0 11 0 11]);
    xticks(0:1:11);
    yticks(0:1:11);
    grid on;
        % Get the current tick labels and replace '11' with ''
    xticklabels = get(gca,'XTickLabel');
    yticklabels = get(gca,'YTickLabel');
    xticklabels{end} = '';
    yticklabels{end} = '';
    
    % Apply the modified tick labels
    set(gca, 'XTickLabel', xticklabels);
    set(gca, 'YTickLabel', yticklabels);
    hold off;
end

% Add a legend outside the subplots
legend({'Action 1', 'Action 2', 'Action 3'}, 'Location', 'outside');


figure;
% Find indices where b1, b2, and b3 are all equal to 3
selected_indices = find(state_space(:,3) == 2 & state_space(:,4) == 3 & state_space(:,5) == 3);

% Extract AoI1, AoI2, and corresponding policy choices for these indices
AoI1_selected = state_space(selected_indices, 1);
AoI2_selected = state_space(selected_indices, 2);
policy_selected = policy(selected_indices);

% Generate the scatter plot with different colors for each action
scatter(AoI1_selected, AoI2_selected, 100, policy_selected, 'filled');
xlabel('AoI1');
ylabel('AoI2');
title('Scatter plot for policy when b1=b2=b3=3');
colorbar;
grid on;
colormap(jet); % Use jet colormap for better differentiation between actions
caxis([1 max(policy)]); % Adjust the color axis to match policy range




