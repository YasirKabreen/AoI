clc;
clear;

% The state variables are [AoI1, AoI2, AoI3, B_level_1, B_level_2, B_level_3, w_ch1, w_ch2, w_ch3].
AoI1 = [1, 2, ];
AoI2 = [1, 2, ];
AoI3 = [1, 2, ];
B_level_1 = [0, 1, 2];
B_level_2 = [0, 1, 2];
B_level_3 = [0, 1, 2];
w_ch1 = [1, 2]; % 1 = Good, 2 = Bad.
w_ch2 = [1, 2]; % 1 = Good, 2 = Bad.
w_ch3 = [1, 2]; % 1 = Good, 2 = Bad.

% Generate all combinations of state variables
StateSpace = [];
for i = 1:length(AoI1)
    for j = 1:length(AoI2)
        for k = 1:length(AoI3)
            for l = 1:length(B_level_1)
                for m = 1:length(B_level_2)
                    for n = 1:length(B_level_3)
                        for o = 1:length(w_ch1)
                            for p = 1:length(w_ch2)
                                for q = 1:length(w_ch3)
                                    StateSpace = [StateSpace; [AoI1(i), AoI2(j), AoI3(k), B_level_1(l), B_level_2(m), B_level_3(n), w_ch1(o), w_ch2(p), w_ch3(q)]];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

x = zeros(length(StateSpace), length(StateSpace));


% when the first Sensor is sending.
n = StateSpace;
n1 = [];
n2 = [];
n3 = [];
x1 = n;
x2 = n;

for i = 1:length(StateSpace)
    x1(i, 7) = 1;
    if (x1(i, 4) > 0) && (x1(i, 7) == 1)
        x1(i, 1) = 1;
        x1(i, 4) = x1(i, 4) - 1;
    elseif (x1(i, 4) > 1) && (x1(i, 7) == 2)
        x1(i, 1) = 1;
        x1(i, 4) = x1(i, 4) - 2;
    else
        x1(i, 1) = min(x1(i, 1) + 1, max(AoI1));
    end

    x1(i, 2) = min(x1(i, 2) + 1, max(AoI2));
    x1(i, 3) = min(x1(i, 3) + 1, max(AoI3));
    n1 = [n1; x1(i, :)];
end

for i = 1:length(StateSpace)
    x2(i, 7) = 2;
    if (x2(i, 4) > 0) && (x2(i, 7) == 1)
        x2(i, 1) = 1;
        x2(i, 4) = x2(i, 4) - 1;
    elseif (x2(i, 4) > 1) && (x2(i, 7) == 2)
        x2(i, 1) = 1;
        x2(i, 4) = x2(i, 4) - 2;
    else
        x2(i, 1) = min(x2(i, 1) + 1, max(AoI1));
    end

    x2(i, 2) = min(x2(i, 2) + 1, max(AoI2));
    x2(i, 3) = min(x2(i, 3) + 1, max(AoI3));
    n2 = [n2; x2(i, :)];
end

% when the second Sensor is sending
m = StateSpace;
m1 = [];
m2 = [];
m3 = [];
x1 = m;
x2 = m;

for i = 1:length(StateSpace)
    x1(i, 8) = 1;
    if (x1(i, 5) > 0) && (x1(i, 8) == 1)
        x1(i, 2) = 1;
        x1(i, 5) = x1(i, 5) - 1;
    elseif (x1(i, 5) > 1) && (x1(i, 8) == 2)
        x1(i, 2) = 1;
        x1(i, 5) = x1(i, 5) - 2;
    else
        x1(i, 2) = min(x1(i, 2) + 1, max(AoI2));
    end

    x1(i, 1) = min(x1(i, 1) + 1, max(AoI1));
    x1(i, 3) = min(x1(i, 3) + 1, max(AoI3));
    m1 = [m1; x1(i, :)];
end

for i = 1:length(StateSpace)
    x2(i, 8) = 2;
    if (x2(i, 5) > 0) && (x2(i, 8) == 1)
        x2(i, 2) = 1;
        x2(i, 5) = x2(i, 5) - 1;
    elseif (x2(i, 5) > 1) && (x2(i, 8) == 2)
        x2(i, 2) = 1;
        x2(i, 5) = x2(i, 5) - 2;
    else
        x2(i, 2) = min(x2(i, 2) + 1, max(AoI2));
    end

    x2(i, 1) = min(x2(i, 1) + 1, max(AoI1));
    x2(i, 3) = min(x2(i, 3) + 1, max(AoI3));
    m2 = [m2; x2(i, :)];
end

% when the third Sensor is sending
o = StateSpace;
o1 = [];
o2 = [];
x1 = o;
x2 = o;

for i = 1:length(StateSpace)
    x1(i, 9) = 1;
    if (x1(i, 6) > 0) && (x1(i, 9) == 1)
        x1(i, 3) = 1;
        x1(i, 6) = x1(i, 6) - 1;
    elseif (x1(i, 6) > 1) && (x1(i, 9) == 2)
        x1(i, 3) = 1;
        x1(i, 6) = x1(i, 6) - 2;
    else
        x1(i, 3) = min(x1(i, 3) + 1, max(AoI3));
    end

    x1(i, 1) = min(x1(i, 1) + 1, max(AoI1));
    x1(i, 2) = min(x1(i, 2) + 1, max(AoI2));
    o1 = [o1; x1(i, :)];
end

for i = 1:length(StateSpace)
    x2(i, 9) = 2;
    if (x2(i, 6) > 0) && (x2(i, 9) == 1)
        x2(i, 3) = 1;
        x2(i, 6) = x2(i, 6) - 1;
    elseif (x2(i, 6) > 1) && (x2(i, 9) == 2)
        x2(i, 3) = 1;
        x2(i, 6) = x2(i, 6) - 2;
    else
        x2(i, 3) = min(x2(i, 3) + 1, max(AoI3));
    end

    x2(i, 1) = min(x2(i, 1) + 1, max(AoI1));
    x2(i, 2) = min(x2(i, 2) + 1, max(AoI2));
    o2 = [o2; x2(i, :)];
end


v1=[];
v2=[];
for i=1 : length(n1)
    for j=1 : length(StateSpace)
        if (n1(i,:)==StateSpace(j,:))
            v1=[v1;j];
        end
        if (n2(i,:)==StateSpace(j,:))
            v2=[v2;j];
        end
    end
end

w1=[];
w2=[];
for i=1 : length(m)
    for j=1 : length(StateSpace)
        if (m1(i,:)==StateSpace(j,:))
            w1=[w1;j];
        end
        if (m2(i,:)==StateSpace(j,:))
            w2=[w2;j];
        end
    end
end

z1=[];
z2=[];
for i=1 : length(o)
    for j=1 : length(StateSpace)
        if (o1(i,:)==StateSpace(j,:))
            z1=[z1;j];
        end
        if (o2(i,:)==StateSpace(j,:))
            z2=[z2;j];
        end
    end
end



p1 = x;
for i = 1:length(v1)
    p1(i, v1(i)) = 0.5;
    p1(i, v2(i)) = 0.5;
end

p2 = x;
for i = 1:length(w1)
    p2(i, w1(i)) = 0.7;
    p2(i, w2(i)) = 0.3;
end

p3 = x;
for i = 1:length(w1)
    p3(i, z1(i)) = 0.6;
    p3(i, z2(i)) = 0.4;
end

P(:, :, 1) = p1;
P(:, :, 2) = p2;
P(:, :, 3) = p3;
R(:, 1) = n(:, 1)';
R(:, 2) = m(:, 2)';
R(:, 3) = o(:, 3)';

mdp_check(P, R);

discount = 0.95;
[V, policy] = mdp_policy_iteration(P, R, discount);
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

% Simulation run using the policy
starting_state_index = 57; % You can change the starting state index here
simulation_time_steps = 100;

current_state = StateSpace(starting_state_index, :);
accumulated_reward = 0;

for time_step = 1:simulation_time_steps
    action = policy(starting_state_index);
    reward = R(starting_state_index, action);
    accumulated_reward = accumulated_reward + reward;
    
    % Get the probabilities for the next state given the current state and action
    prob_next_state = P(starting_state_index, :, action);
    
    % Convert the probabilities to a cumulative distribution function (CDF)
    cdf_next_state = cumsum(prob_next_state);
    
    % Generate a random number to determine the next state
    random_num = rand();
    next_state_index = find(cdf_next_state >= random_num, 1, 'first');
    
    % Update the current state for the next iteration
    current_state = StateSpace(next_state_index, :);
    starting_state_index = next_state_index;
end

fprintf('Total accumulated reward after %d time steps:%f\n', simulation_time_steps, accumulated_reward);

