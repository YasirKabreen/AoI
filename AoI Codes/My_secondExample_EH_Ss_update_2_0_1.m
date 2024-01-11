clc;
clear;

% Define the state space components
AoI1_values = 1:10;
AoI2_values = 1:10;
B = 5; % Maximum battery level
Battery_values = 1:B;

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

for i = 1:num_states
    state = state_space(i, :);
    delta1 = state(1);
    delta2 = state(2);
    b1 = state(3);
    b2 = state(4);
    b3 = state(5);

    % Battery penalty per sensor
    battery_penalty1 = b1 * 0.05;
    battery_penalty2 = b2 * 0.05;
    battery_penalty3 = b3 * 0.05;

    % Cost computation with battery penalty
    Cost1 = (delta1 + delta2)/2 + 1 - 0.5 * p * q * delta1 + battery_penalty1;
    Cost2 = (delta1 + delta2)/2 + 1 - 0.5 * p * q * delta2 + battery_penalty2;
    Cost3 = (delta1 + delta2)/2 + 1 - 0.5 * (1-p) * q * (delta1 + delta2) + battery_penalty3;

    % Rewards
    R1(i) = -Cost1;
    R2(i) = -Cost2;
    R3(i) = -Cost3;

   % New states after considering energy harvesting for each action:

% For successful transmissions:
new_state1_eh = [1, min(delta2+1, 10), min(max(b1-1, 1)+1, B), min(b2+1, B), min(b3+1, B)];
new_state2_eh = [min(delta1+1, 10), 1, min(b1+1, B), min(max(b2-1, 1)+1, B), min(b3+1, B)];
new_state3_eh = [1, 1, min(b1+1, B), min(b2+1, B), min(max(b3-1, 1)+1, B)];

% For unsuccessful transmissions:
new_state_fail1_eh = [min(delta1+1, 10), min(delta2+1, 10), min(max(b1-1, 1)+1, B), min(b2+1, B), min(b3+1, B)];
new_state_fail2_eh = [min(delta1+1, 10), min(delta2+1, 10), min(b1+1, B), min(max(b2-1, 1)+1, B), min(b3+1, B)];
new_state_fail3_eh = [min(delta1+1, 10), min(delta2+1, 10), min(b1+1, B), min(b2+1, B), min(max(b3-1, 1)+1, B)];


% New states after considering no energy harvesting for each action:

% For successful transmissions:
new_state1_no_eh = [1, min(delta2+1, 10), max(b1-1, 1), b2, b3];
new_state2_no_eh = [min(delta1+1, 10), 1, b1, max(b2-1, 1), b3];
new_state3_no_eh = [1, 1, b1, b2, max(b3-1, 1)];

% For unsuccessful transmissions:
new_state_fail1_no_eh = [min(delta1+1, 10), min(delta2+1, 10), max(b1-1, 1), b2, b3];
new_state_fail2_no_eh = [min(delta1+1, 10), min(delta2+1, 10), b1, max(b2-1, 1), b3];
new_state_fail3_no_eh = [min(delta1+1, 10), min(delta2+1, 10), b1, b2, max(b3-1, 1)];


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

% Visualization and plotting remain the same as in the previous snippet...
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
    hold on;
end

% Add labels
xlabel('AoI1');
ylabel('AoI2');
axis([0 21 0 21]);
xticks(0:1:21);
yticks(0:1:21);
legend('Action 1', 'Action 2', 'Action 3');
grid on;
