                                                                                                                                                           
 clc,clear 

% The state variable are [AoI , B_level , w_ch_state].
AoI1=[1 2 3 4 5 6]; 
AoI2=[1 2 3 4 5 6] ;  
B_level_1=[0 1 2 3 4 5];  
B_level_2=[0 1 2 3 4 5]; 
w_ch1 = [1 ,2];  
% 1 = Good , 2 = Bad.

w_ch2 = [1 ,2];  
% 1 = Good , 2 = Bad.
  StateSpace=[];
  for i=1 : length(AoI1)
      for j=1 : length(AoI2)
          for k=1 : length(B_level_1)
             for l=1 : length(B_level_2)
                  for e=1 : length(w_ch1)
                      for ee=1 : length(w_ch2)
                          StateSpace=[StateSpace;[AoI1(i),AoI2(j),B_level_1(k),B_level_2(l),w_ch1(e),w_ch2(ee)]];
                      end
                  end 
             end    
          end
      end
   end
  
StateSpace
x=zeros(length(StateSpace),length(StateSpace));

  % when the first Sensor is sending.
  n=StateSpace;
  n1=[];
  n2=[];
  x1=[];
  x2=[];
  for i =1 : length(StateSpace)
      x1(i,:)=n(i,:);
    x2(i,:)=n(i,:);
      for ii=1:2
         
         if(ii==1)
              x1(i,5)=ii;
      if (x1(i,3)>0)&&(x1(i,5)==1)
         x1(i,1)=1;
         x1(i,3)=x1(i,3)-1;
      elseif(x1(i,3)>1)&&(x1(i,5)==2)
            
          x1(i,1)=1;
          x1(i,3)=x1(i,3)-2;
      else 
      
          x1(i,1)= min(x1(i,1)+1,max(AoI1));
      end
      
     
      
      x1(i,2)= min(x1(i,2)+1,max(AoI2));
           n1=[n1;x1(i,:)];
         else
               x2(i,5)=ii;
      if (x2(i,3)>0)&&(x2(i,5)==1)
         x2(i,1)=1;
         x2(i,3)=x2(i,3)-1;
      elseif(x2(i,3)>1)&&(x2(i,5)==2)
            
          x2(i,1)=1;
          x2(i,3)=x2(i,3)-2;
      else 
      
          x2(i,1)= min(x2(i,1)+1,max(AoI1));
      end
      
     
      
      x2(i,2)= min(x2(i,2)+1,max(AoI2));
             n2=[n2;x2(i,:)];
         end
      end
  end
  
  
  
  % when the second Sensor is sending
 m=StateSpace;
 m1=[];
 m2=[];
 x1=[];
  x2=[];
  for i =1 : length(StateSpace)
      x1(i,:)=m(i,:);
    x2(i,:)=m(i,:);
      for ii=1:2
         
         if(ii==1)
              x1(i,6)=ii;
      if (x1(i,4)>0)&&(x1(i,6)==1)
         x1(i,2)=1;
         x1(i,4)=x1(i,4)-1;
      elseif(x1(i,4)>1)&&(x1(i,6)==2)
            
          x1(i,2)=1;
          x1(i,4)=x1(i,4)-2;
      else 
      
          x1(i,2)= min(x1(i,2)+1,max(AoI2));
      end
      
     
      
      x1(i,1)= min(x1(i,1)+1,max(AoI1));
           m1=[m1;x1(i,:)];
         else
               x2(i,6)=ii;
      if (x2(i,4)>0)&&(x2(i,6)==1)
         x2(i,2)=1;
         x2(i,4)=x2(i,4)-1;
      elseif(x2(i,4)>1)&&(x2(i,6)==2)
            
          x2(i,2)=1;
          x2(i,4)=x2(i,4)-2;
      else 
      
          x2(i,2)= min(x2(i,2)+1,max(AoI2));
      end
      
     
      
      x2(i,1)= min(x2(i,1)+1,max(AoI1));
             m2=[m2;x2(i,:)];
         end
      end
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

p1=x;
for i = 1 : length(v1)
p1(i,v1(i))=0.5;
p1(i,v2(i))=0.5;
end

p2=x;
for i = 1 : length(w1)
p2(i,w1(i))=0.7;
p2(i,w2(i))=0.3;
end


 MDPtoolbox_path = pwd;
addpath(MDPtoolbox_path)
P(:,:,1) = p1;
P(:,:,2) = p2;
 R(:,1) =3-n(:,1)'-n(:,3)'  ;
R(:,2) =3-m(:,2)'-n(:,4)' ;

mdp_check(P, R)

discount = 0.95;
[V, policy] = mdp_policy_iteration(P, R, discount)


% Simulation run using the policy
starting_state_index = 141; % Modify as needed
simulation_time_steps = 100;
accumulated_reward = 0;

current_state_index = starting_state_index;

% Lists to store AoI values for each sensor
aoi_sensor_1 = zeros(1, simulation_time_steps);
aoi_sensor_2 = zeros(1, simulation_time_steps);


for time_step = 1:simulation_time_steps
    action = policy(current_state_index);
    reward = R(current_state_index, action);
    accumulated_reward = accumulated_reward + reward;

    % Get the probabilities for the next state given the current state and action
    prob_next_state = P(current_state_index, :, action);

    % Convert the probabilities to a cumulative distribution function (CDF)
    cdf_next_state = cumsum(prob_next_state);

    % Generate a random number to determine the next state
    random_num = rand();
    next_state_index = find(cdf_next_state >= random_num, 1, 'first');

    % Update the current state index for the next iteration
    current_state_index = next_state_index;

    % Check if current_state_index is still within valid range
    if current_state_index <= size(StateSpace, 1)
        % Calculate AoI for each sensor
        aoi_sensor_1(time_step) = StateSpace(current_state_index, 1);
        aoi_sensor_2(time_step) = StateSpace(current_state_index, 2);
    else
        % Handle the case where current_state_index is out of range
        % You might want to break the loop or take some other action
    end
end


figure;
subplot(2, 1, 1);
plot(policy);
xlabel('State Index');
ylabel('Action');
title('Optimal Policy for AoI Minimization');
axis([0 length(StateSpace)+1 0 4]);
xticks(1:length(StateSpace));
yticks(0:4);
yticklabels({'No Action', 'Sensor 1', 'Sensor 2', 'Sensor 3'});
grid on;

% Plotting AoI values for each sensor
subplot(2, 1, 2);
plot(1:simulation_time_steps, aoi_sensor_1, 'b', 1:simulation_time_steps, aoi_sensor_2, 'g');
xlabel('Time Step');
ylabel('AoI');
title('AoI Evolution for Each Sensor');
legend('Sensor 1', 'Sensor 2');
axis([0 simulation_time_steps 0 7]);
grid on;
grid minor;

