clc;
clear;

% Define the state space components
AoI1_values = (1:20);
AoI2_values = (1:20);
Ss1=[0 1];
Ss2=[0 1];



% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values) * numel(Ss1)* numel(Ss2);

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 4); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
        for Ss1 = Ss1
            for Ss2 = Ss2
        
                        state_space(idx, :) = [aoi1, aoi2,Ss1,Ss2];
                        idx = idx + 1;
            end
        end
    end
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);



% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);



discount = 0.95;

% Loop through all states to calculate transitions and rewards

p = input('enter the probability p'); % Probability of sending
q = input('enter the probability q'); % Probability of source state





ok the logic are :

q is probability of Source state 
p is prbability of ttransmission

 if Ss ==0 
with probability q    { A=min(A+1,max(A)); , Ss=Ss';}
with probability 1-q  { A=min(A+1,max(A)); , Ss=SS;}


if Ss==1

with probability q * p  {A=1; , Ss=Ss' ;}

with probability 1-q * p  { A=1; , Ss=SS;}

with probability q * 1-p  {A=min(A+1,max(A)); , Ss=Ss';}

with probability 1-q * 1-p  {A=min(A+1,max(A)); , Ss=Ss;}

