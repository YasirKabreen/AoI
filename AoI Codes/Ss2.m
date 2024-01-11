clc,clear
% Define the state space components
AoI1_values = [1, 2,];
AoI2_values = [1, 2,];
LocationAGV1_values = [1, 2, 3, 4];
LocationAGV2_values = [1, 2, 3, 4];
Source1State_values = [1];
Source2State_values = [1];
% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * ...
numel(LocationAGV1_values) * numel(LocationAGV2_values) * ...
numel(Source1State_values) * numel(Source2State_values);
% Create a matrix to store all combinations of states
% Each row represents a combination of state values
state_space = zeros(num_states, 6); % 6 columns for state components
% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
for aoi2 = AoI2_values
for locAGV1 = LocationAGV1_values
for locAGV2 = LocationAGV2_values
for source1 = Source1State_values
for source2 = Source2State_values
state_space(idx, :) = [aoi1, aoi2, locAGV1, locAGV2, source1, source2];
idx = idx + 1;
end
end
end
end
end
end
% Now, state_space is a matrix containing all combinations of state values
state_space;
%sensor ont to send
n= state_space;
for i=1 :length(state_space)
n(i,1)= 1;
if (n(i,3)==4)
n(i,3)=1;
else
n(i,3)=n(i,3)+1;
end
n(i,2)=min(n(i,2)+1 ,max(AoI2_values));
if (n(i,4)==4)
n(i,4)=1;
else
n(i,4)=n(i,4)+1;
end
end
n
%sensor two to send
m=state_space;
for i=1 :length(state_space)
m(i,2)= 1;
if (m(i,4)==4)
m(i,4)=1;
else
m(i,4)=m(i,4)+1;
end
m(i,1)=min(m(i,1)+1 ,max(AoI1_values));
if (m(i,3)==4)
m(i,3)=1;
else
m(i,3)=m(i,3)+1;
end
end
m
v=[];
for i=1 : length(n)
    for j=1 : length(state_space)
        if (n(i,:)==state_space(j,:))
            v=[v;j];
        end
        
    end
end
v

w=[];
for i=1 : length(m)
    for j=1 : length(state_space)
        if (m(i,:)==state_space(j,:))
            w=[w;j];
        end
        
    end
end
w
x = zeros(length(state_space), length(state_space));
p1 = x;
for i = 1:length(v)
    p1(i, v(i)) = 1;
  
end
p1


p2 = x;
for i = 1:length(w)
    p1(i, w(i)) = 1;
  
end

P(:, :, 1) = p1;
P(:, :, 2) = p2;

R(:, 1) = 1./n(:, 1)';
R(:, 2) =1./ m(:, 2)';


mdp_check(P, R);

discount = 0.95;
[V, policy] = mdp_policy_iteration(P, R, discount)










