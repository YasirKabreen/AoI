% Function to simulate a given policy
function [sim_states, sim_actions, sim_rewards] = simulate_policy(P, R, state_space, initial_state, num_steps, policy)
    sim_states = zeros(num_steps, 5);
    sim_actions = zeros(num_steps, 1);
    sim_rewards = zeros(num_steps, 1);
    
    num_states = size(state_space, 1);

    sim_states(1, :) = initial_state;

    for step = 1:num_steps
        current_state = sim_states(step, :);
        state_idx = find(ismember(state_space, current_state, 'rows'));

        if isempty(policy)
            action = randi(4);  % Choose action uniformly at random
        else
            action = policy(state_idx);
        end

        sim_actions(step) = action;
        next_state_idx = randsample(num_states, 1, true, P(state_idx,:,action));
        next_state = state_space(next_state_idx, :);

        sim_states(step+1, :) = next_state;
        sim_rewards(step) = R(state_idx, action);
    end
end