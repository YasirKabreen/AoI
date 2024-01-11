clc,clear
R=[.1 .8 0 .1 0 0 0 0; 0 .2 .8 0 0 0 0 0; 0 0 .9 0 .1 0 0 0; .1 0 0 .8 0 .1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 .1 0 .1 .8 0; 0 0 0 0 0 0 .2 .8; 0 0 0 0 0 0 0 0];
L=[.9 0 0 .1 0 0 0 0; .8 .2 0 0 0 0 0 0; 0 .8 .1 0 .1 0 0 0; .1 0 0 .8 0 .1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 .1 0 .9 0 0; 0 0 0 0 0 .2 .8 0; 0 0 0 0 0 0 0 0];
UP=[.1 .1 0 .8 0 0 0 0; .1 .8 .1 0 0 0 0 0; 0 .1 .1 0 .8 0 0 0; 0 0 0 .2 0 .8 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 .9 .1 0; 0 0 0 0 0 .1 .8 .1; 0 0 0 0 0 0 0 0];
DW=[.9 .1 0 0 0 0 0 0; .1 .8 .1 0 0 0 0 0; 0 .1 .9 0 0 0 0 0; .8 0 0 .2 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 .8 0 .1 .1 0; 0 0 0 0 0 .1 .8 .1; 0 0 0 0 0 0 0 0];
reward=[0 0 0 0 -1 0 0 1]';
V=[0 0 0 0 0 0 0 0]';
Vnext=[0 0 0 0 0 0 0 0]';
g=0.9;

for i=1:100
VnextL=reward+g*L*V;
VnextR=reward+g*R*V;
VnextUP=reward+g*UP*V;
VnextDW=reward+g*DW*V;

for j=1:length(V)
   Vnext(j)= max(max([VnextL(j),VnextR(j),VnextUP(j),VnextDW(j)])); 
   
if VnextL(j) == Vnext(j)
      policy(j)='L';
elseif VnextR(j) == Vnext(j)
      policy(j)='R';
elseif VnextUP(j) == Vnext(j)
      policy(j)='U';
else
  policy(j)='D';

end

end

V=Vnext;

end
V
policy'



  
 
 
 
 
 
 
 
 
 