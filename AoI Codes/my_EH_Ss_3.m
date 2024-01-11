clc;
clear;

% Define the state space components
AoI1_values = (1:4);
AoI2_values = (1:4);

B_level1=(0:4);
B_level2=(0:4);
B_level3=(0:4);

% Calculate the number of states
num_states = numel(AoI1_values) * numel(AoI2_values)*   numel(B_level1) * numel(B_level2)*   numel(B_level3)  ;

% Create a matrix to store all combinations of states 
state_space = zeros(num_states, 5); 

% Generate all combinations of states
idx = 1;
for aoi1 = AoI1_values
    for aoi2 = AoI2_values
         for b1 = B_level1
            for b2 = B_level2
                 for b3 = B_level3
                        state_space(idx, :) = [aoi1, aoi2,b1,b2,b3];
                        idx = idx + 1;
                 end
            end
         end
    end
end

% Initialize transition probability matrices
P1 = zeros(num_states, num_states);
P2 = zeros(num_states, num_states);
P3 = zeros(num_states, num_states);
P4 = zeros(num_states, num_states);


% Initialize reward matrix
R1 = zeros(num_states, 1);
R2 = zeros(num_states, 1);
R3 = zeros(num_states, 1);
R4 = zeros(num_states, 1);


discount = 0.95;

% Loop through all states to calculate transitions and rewards

e= input('enter the probability p'); % Probability of EH
q = input('enter the probability q'); % Probability of source state
for i = 1:num_states
    state = state_space(i, :);
    
  % when send Sensor one transitions and rewards
    next_state11 = state;
    next_state12 = state;
    next_state13 = state;
    next_state14 = state;
  %with probability e*q
  if state(3)>=1 
        next_state11(1) = 1;
        next_state11(2)=min(next_state11(2)+1,max(AoI2_values));
        next_state11(3)=next_state11(3);
        next_state11(4)=min(next_state11(4)+1,max(B_level2));
        next_state11(5)=min(next_state11(5)+1,max(B_level3));
   
       
  else
         next_state11(1) = min(next_state11(1) + 1, max(AoI1_values));
         next_state11(2) = min(next_state11(2) + 1, max(AoI2_values));
         next_state11(3) = min(next_state11(3)+1,max(B_level1));
         next_state11(4) = min(next_state11(4)+1,max(B_level2));
         next_state11(5) = min(next_state11(5)+1,max(B_level3));
  end
  
    %with probability 1-e*q
      if state(3)>=1 
        next_state12(1) = 1;
        next_state12(2)=min(next_state12(2)+1,max(AoI2_values));
        next_state12(3)=next_state12(3)-1;
        next_state12(4)=next_state12(4);
        next_state12(5)=next_state12(5);
   
       
        else
         next_state12(1) = min(next_state12(1) + 1, max(AoI1_values));
         next_state12(2) = min(next_state12(2) + 1, max(AoI2_values));
         next_state12(3) = next_state12(3);
         next_state12(4) = next_state12(4);
         next_state12(5) = next_state12(5);
       end
    
      %with probability e*1-q
        next_state13(1) = min(next_state13(1) + 1, max(AoI1_values));
        next_state13(2) = min(next_state13(2) + 1, max(AoI2_values));
        next_state13(3)=min(next_state13(3)+1,max(B_level1));
        next_state13(4)=min(next_state13(4)+1,max(B_level2));
        next_state13(5)=min(next_state13(5)+1,max(B_level3));
    
      %with probability 1-e*1-q
        next_state14(1) = min(next_state14(1) + 1, max(AoI1_values));
        next_state14(2) = min(next_state14(2) + 1, max(AoI2_values));
        
        
        
  % when send Sensor two transitions and rewards
    next_state21 = state;
    next_state22 = state;
    next_state23 = state;
    next_state24 = state;
  %with probability e*q
  if state(4)>=1 
        next_state21(2) = 1;
        next_state21(1)=min(next_state21(1)+1,max(AoI1_values));
        next_state21(4)=next_state21(4);
        next_state21(3)=min(next_state21(3)+1,max(B_level1));
        next_state21(5)=min(next_state21(5)+1,max(B_level3));
   
       
  else
         next_state21(1) = min(next_state21(1) + 1, max(AoI1_values));
         next_state21(2) = min(next_state21(2) + 1, max(AoI2_values));
         next_state21(3) = min(next_state21(3)+1,max(B_level1));
         next_state21(4) = min(next_state21(4)+1,max(B_level2));
         next_state21(5) = min(next_state21(5)+1,max(B_level3));
  end
  
    %with probability 1-e*q
      if state(4)>=1 
        next_state22(2) = 1;
        next_state22(1)=min(next_state22(1)+1,max(AoI1_values));
        next_state22(4)=next_state22(4)-1;
        next_state22(3)=next_state22(3);
        next_state22(5)=next_state22(5);
   
       
        else
         next_state22(1) = min(next_state22(1) + 1, max(AoI1_values));
         next_state22(2) = min(next_state22(2) + 1, max(AoI2_values));
         next_state22(3) = next_state22(3);
         next_state22(4) = next_state22(4);
         next_state22(5) = next_state22(5);
       end
    
      %with probability e*1-q
        next_state23(1) = min(next_state23(1) + 1, max(AoI1_values));
        next_state23(2) = min(next_state23(2) + 1, max(AoI2_values));
        next_state23(3)=min(next_state23(3)+1,max(B_level1));
        next_state23(4)=min(next_state23(4)+1,max(B_level2));
        next_state23(5)=min(next_state23(5)+1,max(B_level3));
    
      %with probability 1-e*1-q
        next_state24(1) = min(next_state24(1) + 1, max(AoI1_values));
        next_state24(2) = min(next_state24(2) + 1, max(AoI2_values));
     
  % when send Sensor three transitions and rewards
    next_state31 = state;
    next_state32 = state;
    next_state33 = state;
    next_state34 = state;
    

  %with probability e*1-q
  if state(5)>=1 
        next_state31(1) = 1;
        next_state31(2) = 1;
        next_state31(3)=min(next_state31(3)+1,max(B_level1));
        next_state31(4)=min(next_state31(4)+1,max(B_level2));
        next_state31(5)=next_state31(5);
   
       
  else
         next_state31(1) = min(next_state31(1) + 1, max(AoI1_values));
         next_state31(2) = min(next_state31(2) + 1, max(AoI2_values));
         next_state31(3) = min(next_state31(3)+1,max(B_level1));
         next_state31(4) = min(next_state31(4)+1,max(B_level2));
         next_state31(5) = min(next_state31(5)+1,max(B_level3));
  end
  
   %with probability 1-e*1-q
      if state(5)>=1 
        next_state32(1) = 1;
        next_state32(2) = 1;
        next_state32(3)=next_state32(3);
        next_state32(4)=next_state32(4);
        next_state32(5)=next_state32(5)-1;
   
       
        else
         next_state32(1) = min(next_state32(1) + 1, max(AoI1_values));
         next_state32(2) = min(next_state32(2) + 1, max(AoI2_values));
         next_state32(3) = next_state32(3);
         next_state32(4) = next_state32(4);
         next_state32(5) = next_state32(5);
       end
    
       %with probability e*q
        next_state33(1) = min(next_state33(1) + 1, max(AoI1_values));
        next_state33(2) = min(next_state33(2) + 1, max(AoI2_values));
        next_state33(3)=min(next_state33(3)+1,max(B_level1));
        next_state33(4)=min(next_state33(4)+1,max(B_level2));
        next_state33(5)=min(next_state33(5)+1,max(B_level3));
    
      
       %with probability 1-e*q
        next_state34(1) = min(next_state34(1) + 1, max(AoI1_values));
        next_state34(2) = min(next_state34(2) + 1, max(AoI2_values));
        
  % when not send  transitions and rewards
    next_state41 = state;
    next_state42 = state;
   
       %with probability e
         next_state41(1) = min(next_state41(1) + 1, max(AoI1_values));
         next_state41(2) = min(next_state41(2) + 1, max(AoI2_values));
         next_state41(3) = min(next_state41(3)+1,max(B_level1));
         next_state41(4) = min(next_state41(4)+1,max(B_level2));
         next_state41(5) = min(next_state41(5)+1,max(B_level3));
 
 
       %with probability 1-e
        next_state42(1) = min(next_state42(1) + 1, max(AoI1_values));
        next_state42(2) = min(next_state42(2) + 1, max(AoI2_values));      
    
   
    % Find indices of next states in state_space
    next_idx11 = find(ismember(state_space, next_state11, 'rows'));
       next_idx12 = find(ismember(state_space, next_state12, 'rows'));
          next_idx13 = find(ismember(state_space, next_state13, 'rows'));
             next_idx14 = find(ismember(state_space, next_state14, 'rows'));
             
    next_idx21 = find(ismember(state_space, next_state21, 'rows'));
       next_idx22 = find(ismember(state_space, next_state22, 'rows'));
          next_idx23 = find(ismember(state_space, next_state23, 'rows'));
             next_idx24 = find(ismember(state_space, next_state24, 'rows'));
             
    next_idx31 = find(ismember(state_space, next_state31, 'rows'));
       next_idx32 = find(ismember(state_space, next_state32, 'rows'));
          next_idx33 = find(ismember(state_space, next_state33, 'rows'));
             next_idx34 = find(ismember(state_space, next_state34, 'rows'));
             
    next_idx41 = find(ismember(state_space, next_state41, 'rows'));
       next_idx42 = find(ismember(state_space, next_state42, 'rows'));
          
    % Update transition probability matrices
    P1(i, next_idx11) = e*q;
    P1(i, next_idx12) = (1-e)*q;
    P1(i, next_idx13) = e*(1-q);
    P1(i, next_idx14) = (1-e)*(1-q);
    
    P2(i, next_idx21) = e*q;
    P2(i, next_idx22) = (1-e)*q;
    P2(i, next_idx23) = e*(1-q);
    P2(i, next_idx24) = (1-e)*(1-q);
        
    P3(i, next_idx31) = e*(1-q);
    P3(i, next_idx32) = (1-e)*(1-q);
    P3(i, next_idx33) = e*q;
    P3(i, next_idx34) = (1-e)*q;
   
    P4(i, next_idx41) = e ;
    P4(i, next_idx42) = (1-e) ;
 

   
    
     % Update reward matrices
     aa11=(next_state11(1)+next_state12(1)+next_state13(1)+next_state14(1))/4;
     aa12=(next_state11(2)+next_state12(2)+next_state13(2)+next_state14(2))/4;
 
     aa21=(next_state21(1)+next_state22(1)+next_state23(1)+next_state24(1))/4;
     aa22=(next_state21(2)+next_state22(2)+next_state23(2)+next_state24(2))/4;
     
     aa31=(next_state31(1)+next_state32(1)+next_state33(1)+next_state34(1))/4;
     aa32=(next_state31(2)+next_state32(2)+next_state33(2)+next_state34(2))/4;
     
     aa41=(next_state41(1)+next_state42(1))/2;
     aa42=(next_state41(2)+next_state42(2))/2;
     
     bb13=(next_state11(3)+next_state12(3)+next_state13(3)+next_state14(3))/4;
     bb14=(next_state11(4)+next_state12(4)+next_state13(4)+next_state14(4))/4;
     bb15=(next_state11(5)+next_state12(5)+next_state13(5)+next_state14(5))/4;
     
     bb23=(next_state21(3)+next_state22(3)+next_state23(3)+next_state24(3))/4;
     bb24=(next_state21(4)+next_state22(4)+next_state23(4)+next_state24(4))/4;
     bb25=(next_state21(5)+next_state22(5)+next_state23(5)+next_state24(5))/4;
      
     bb33=(next_state31(3)+next_state32(3)+next_state33(3)+next_state34(3))/4;
     bb34=(next_state31(4)+next_state32(4)+next_state33(4)+next_state34(4))/4;
     bb35=(next_state31(5)+next_state32(5)+next_state33(5)+next_state34(5))/4;
     
     bb43=(next_state41(3)+next_state42(3))/2;
     bb44=(next_state41(4)+next_state42(4))/2;
     bb45=(next_state41(5)+next_state42(5))/2;
     
     
     
     
    bb2=(next_state21(2)+next_state22(2))/2;
    R1(i) = (bb13+bb14+bb15) /(aa11+aa12);
    R2(i) = (bb23+bb24+bb25)/ (aa21+aa22);
    R3(i) = (bb33+bb34+bb35)/ (aa31+aa32);
    R4(i) = (bb43+bb44+bb45)/ (aa41+aa42);
    
end

% Combine transition probability matrices for all sensors
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;
P(:,:,4) = P4;


% Combine reward matrices for all sensors
R = [R1, R2,R3,R4];

% Perform value iteration to find the optimal value function and policy
[V, policy] = mdp_policy_iteration(P, R, discount);

% Plotting the Policy
figure;
plot(policy);
xlabel('State Index');
ylabel('Action');
title('Optimal Policy for AoI Minimization');
axis([0 length(state_space) + 1 0 4]); % Adjust the axis based on the number of sensors
yticks(1:5);
yticklabels({ 'Sensor 1', 'Sensor 2', 'Sensor 3','No Action',});
grid on;





