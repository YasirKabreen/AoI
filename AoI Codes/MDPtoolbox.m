clc,clear
MDPtoolbox_path = pwd;
addpath(MDPtoolbox_path)
P(:,:,1) = [0 0 0 0 0 0; 0 0 0 0 0 0; 0.3 0.7 0 0 0 0; 0 0 0 0 0 0; 0 0 0.3 0.7 0 0; 0.3 0.7 0 0 0 0];
P(:,:,2) = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
%P(:,:,3) = [0 0 0 0 0 0; 0.1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0.1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0.1 0];
R(:,1) = [1 1 1 1 0 0];
R(:,2) = [0 0 0 0 0 0];
mdp_check(P, R)

discount = 0.95;
[V, policy] = mdp_policy_iteration(P, R, discount)