% ... (Your existing code)

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

% Find the indices of states where B_level2 is equal to 6
target_B_level2 = 6;
target_indices = find(state_space(:, 3) == target_B_level2);

% Extract the corresponding states and policy values for the target indices
target_states = state_space(target_indices, :);
target_policy = policy(target_indices);

% Display the extracted states and policy values
disp('States where B_level2 is equal to 6:');
disp(target_states);
disp('Policy values for these states:');
disp(target_policy);

% ... (Your existing code for plotting)
