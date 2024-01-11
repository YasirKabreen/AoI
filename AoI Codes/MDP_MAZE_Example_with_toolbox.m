clc,clear

discount=0.9;
P(:,:,1)=[.1 .8 0 .1 0 0 0 0; 0 .2 .8 0 0 0 0 0; 0 0 .9 0 .1 0 0 0; .1 0 0 .8 0 .1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 .1 0 .1 .8 0; 0 0 0 0 0 0 .2 .8; 0 0 0 0 0 0 0 0];
P(:,:,2)=[.9 0 0 .1 0 0 0 0; .8 .2 0 0 0 0 0 0; 0 .8 .1 0 .1 0 0 0; .1 0 0 .8 0 .1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 .1 0 .9 0 0; 0 0 0 0 0 .2 .8 0; 0 0 0 0 0 0 0 0];
P(:,:,3)=[.1 .1 0 .8 0 0 0 0; .1 .8 .1 0 0 0 0 0; 0 .1 .1 0 .8 0 0 0; 0 0 0 .2 0 .8 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 .9 .1 0; 0 0 0 0 0 .1 .8 .1; 0 0 0 0 0 0 0 0];
P(:,:,4)=[.9 .1 0 0 0 0 0 0; .1 .8 .1 0 0 0 0 0; 0 .1 .9 0 0 0 0 0; .8 0 0 .2 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 .8 0 .1 .1 0; 0 0 0 0 0 .1 .8 .1; 0 0 0 0 0 0 0 0];
reward=[0 0 0 0 -1 0 0 1]';





    R1=reward;
    R2=reward;
    R3=reward;
    R4=reward;
    R = [R1, R2,R3,R4];

 Q = mdp_Q_learning(P, R, discount, 100000);
 
% Initialize the policy as an array of zeros with the same number of rows as the Q-table
policy = zeros(size(Q, 1), 1);

% Extract the optimal action for each state based on the Q-values
for state = 1:size(Q, 1)
    [~, policy(state)] = max(Q(state, :));
end

% Display the policy
policy






