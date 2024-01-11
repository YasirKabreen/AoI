% Q-learning parameters
Q = zeros(num_states, 2); % Initialize Q-values
gamma = 0.95; % discount factor

% Adaptive parameters
delta_d = 10^-7;
slots = 10^6; % for this example, assuming total slots is 10^8 for Q-learning
epsilon = zeros(1, slots);
alpha = zeros(1, slots);

for t = 1:slots
    % Update epsilon (exploration rate)
    epsilon(t) = 0.02 + 0.98 * exp(-delta_d * t);
    
    % Update alpha (learning rate)
    if t <= 10^7
        alpha(t) = 0.5;
    else
        alpha(t) = 0.01;
    end
end

% Q-learning loop
for t = 1:slots
    % Select current state randomly
    current_state_idx = randi([1 num_states]);
    
    % Exploration-exploitation decision
    if rand() < epsilon(t)
        % Exploration: choose a random action
        action = randi([1 2]);
    else
        % Exploitation: choose the action with maximum Q-value
        [~, action] = max(Q(current_state_idx, :));
    end
    
    % Transition to the next state based on chosen action
    transition_probabilities = P(current_state_idx,:,action);
    next_state_idx = randsample(1:num_states, 1, true, transition_probabilities);
    
    % Q-value update
    reward = R(current_state_idx, action);
    maxQ_next = max(Q(next_state_idx, :));
    Q(current_state_idx, action) = Q(current_state_idx, action) + alpha(t) * (reward + gamma * maxQ_next - Q(current_state_idx, action));
end

% Extract policy from Q-values
[~, policy] = max(Q, [], 2);
