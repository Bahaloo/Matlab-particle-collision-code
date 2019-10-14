clc
clear all
% CD to the input file location
cd 'C:\Users\bahal\Desktop'
if exist('Bodydata') ~=2
    msgID = 'myComponent:inputError';
    msgtext = 'Input does not exist or does not have the expected format.';
    
    ME = MException(msgID,msgtext)
end
% load the input file
Bfile =  fopen('BodyData.txt');
% Bfile = uigetdir('BodyData.txt');
[X0,Y0,Vx,Vy,R] =InputData(Bfile);

%incremental time steps to update the positions and check for collisons
dt = 1e-4;
t = 0; %start time
T = 2;  % Simulation time

j = 1;  %counter
X(:,j) = X0;
Y(:,j) = Y0;
while t <= T
    j = j+1;
    % updating the positions by given velocities
    X(:,j) = X(:,j-1) + Vx*dt;
    Y(:,j) = Y(:,j-1) + Vy*dt;
    
    %checking for collision for particle 1 to N by checking the distances
    % between the spheres, if Collison: calls the collison outcomes from
    % CollisionCalculator.m 
    
    for m = 1:size(X,1)
        for n =  setdiff(m:size(X,1), m)
            if sqrt((X(n,j)-X(m,j))^2+((Y(n,j)-Y(m,j))^2)) <= R(n) + R(m)
                
            end
        end
    end
    t = t+dt;
end

