% This program is written by Hassan Bahaloo Horeh to simulate collision of
% multiple particles, Date: 11/28/2018
% copyright by Hassan Bahaloo (bahaloohoreh.h@husky.neu.edu) 
%****************************** USAGE NOTES *******************************
%%% Change the file address in line 30
%%% change the file name in line 32
%% Run the program, Enjoy!
clc;
clear all
format long
t_end = 8800;                % simulation whole time.
e = 1;                      % coefficient of resilience.  0<e <=1;  not 0!
%% ========================================================================
%   reads the input file and pre processes it
%% ========================================================================
%Usage:
%create the input file in the *.txt format and save it on your computer 
%(for example: MyFile.txt located at C:\Uders\...\Desktop.
%The information in the text file must be as follows:
%First line, a heard of 6 columns separated by comma:
%nr, x0, y0, vx, vy, r
%next lines the numerical values corresponding to each 
%particle as guided in the header line,
%for example, this can be a valid text file:
%nr, x0, y0, vx, vy, r
%1, 0.2, 2, 0, -4, 0.1
%2, 2.3, -1, 1, 3, 0.2
%3, 0, 5, 0, 0, 0.1
%*** note: the particles should not have interreference in advance.
%% ========================================================================
% cd 'C:\Users\bahal\Desktop'           % CD to the location of input file.
% %load the input file
% Bfile =  fopen('particles.txt');
% %calls InputData function to read the data
% [X0,Y0,Vx,Vy,radii] =InputData(Bfile);

%% =============== Create axis for graphical presentations ================

xlim([-25 15]);ylim([-25 15]);     % fix the axis limits.
axis square                        % set the axis aspect ratio to 1:1.
X_lim = xlim;                      % get the axis limits and store in X_lim
Y_lim = ylim;
NumOfParticles = 10;

% calculates the center coordinate of axis
X_center = (X_lim(2)+X_lim(1))/2;
Y_center = (Y_lim(2)+Y_lim(1))/2;

%%  inputs random particle data from a function
[X0, Y0, Vx, Vy, radii] = RandParticles(NumOfParticles, ...
    X_lim,Y_lim,X_center ,Y_center)
%% created the walls as VERY large spheres and adds them to particle list
R_wall = 1e10;
X0 = [X0; R_wall + X_lim(2);        ...
    X_center; -R_wall + X_lim(1); X_center];
Y0 = [Y0;          Y_center; R_wall+Y_lim(2);       ...
              Y_center; -R_wall + Y_lim(1)];
Vx = [Vx; zeros(4,1)];
Vy = [Vy; zeros(4,1)];
radii = [radii; R_wall*ones(4,1)];

NumOfParticles  = length(X0);            % obtaines the number of particles
X = X0;                                  % X location of particles
Y = Y0;                                  % Y location of particles
% checking for initial interference
for m = 1:NumOfParticles -4
    for n =  setdiff(m:NumOfParticles -4, m)
        if abs(round(sqrt((X0(m)-X0(n))^2+(Y0(m)-Y0(n))^2)...
                - (radii(m) + radii(n)),3))==0
            disp('Error: particles are already indented to each other');
            centers = [X0 Y0];
            viscircles(centers,radii,'LineWidth',0.5);% Display the circles.
            return
        end
    end
end
%==========================================================================
% this section calculates the expected time-to-hit (t_hit) between any 2
% mutualmparticles, and finds the min of them (t_hit_mit) to use in a
% priority queue. It calculates the t_hits and stores them in H matrix,
H = inf*ones(NumOfParticles ,NumOfParticles );
for m = 1:NumOfParticles -1
    for n =  setdiff(m:NumOfParticles , m)
        dx =  X(m,:)- X(n,:);
        dy =  Y(m,:)- Y(n,:);
        dVx = Vx(m) - Vx(n);
        dVy = Vy(m) - Vy(n);
        
        dR = [dx dy];
        dV = [dVx dVy];
        
        if dot(dV,dR) < 0
            droot = (dot(dV,dR))^2-dot(dV,dV)*(dot(dR,dR)...
                -(radii(m)+radii(n))^2);
            
            if droot <0; H(m,n) = inf;
            else
                H(m,n) = (-dot(dV,dR)-sqrt(droot))/dot(dV,dV);
            end
            if dot(dV,dV) <= eps; H(m,n) = inf;end
        end
    end
end
% finds the minimum time-to-hit which is the sooest collision time!
t_hitParticle_min = min(min(H));
% finds the lable of particles associating in the collision
[m0,n0] = find(H==t_hitParticle_min);
% specifies a measure time as the hit time which will be updated in each
% collision
t_0 = t_hitParticle_min;
% number of time divisions between each collision as 2^(ndivision)
ndivision = 1;
% initial time step
dt(1) = e*t_hitParticle_min/2^ndivision;
% checks that dt is < infinity
if dt(1) == inf
    dt(1) = 0.5*min(radii)/(max(max(abs(Vx),max(abs(Vy)))))
end
%%
NumOfCollision = 0;                       % will count collision numbers
t_cont=zeros(NumOfParticles ,1);
t = 0;                                     % initial time
k=1;
% constructs a loop to run the simulation to the end time
while t<= t_end
    cla                                    % clears the axes.
        title(['t= ' num2str(t),num2str(t_end),'  Event Driven Sim.'])
        xlabel(['Number of Collisions: ',num2str(NumOfCollision)])
%     title(['Event Driven Simulation'])
    
    % plots a rectangle for just graphical purposes
    Rec = rectangle('Position',[X_lim(1) Y_lim(1) X_lim(2)-X_lim(1) ...
        Y_lim(2)-Y_lim(1)],'LineWidth',0.5);

    % updates the center points of the particles
    centers = [X0 Y0]+[Vx Vy].*(t*ones(NumOfParticles ,1)-t_cont);
    X = centers(:,1);
    Y = centers(:,2);
    % displays the circles.
 
     VISCIR = viscircles(centers,radii,'LineWidth',0.5,'Color','b');
     xticklabels({})
     yticklabels({})
     %======================================================================
    %======================================================================
    % a counter for simultaneous collisions of multiple particles
    NumOfSimCont = 1;
    if abs(round((t- t_0),3)) < eps     % if it is a collision time!
        NumOfCollision = NumOfCollision +1;
        % SimentansCollicionCounter is the counter on simaltaneous 
        %collisions, eg. if 1,2 and 3,7
        % are in collision at the same t, SimentansCollicionCounter = 2;
        for SimentansCollicionCounter = 1:NumOfSimCont

            m0 = m0(SimentansCollicionCounter);
            n0 = n0(SimentansCollicionCounter);
            % mass ratio of the colliding particles assuming spherical
            % shapes and same densities
            alpha = (radii(n0)./radii(m0)).^3;
            
            %coputes the normal vector of the collision and normalizes it
            normal = [(X(m0)-X(n0)),(Y(m0) - Y(n0))];
            normal = normal/norm(normal);
            % normal vector along X direction ex:
            ex = [1,0];
            Vm0 = [Vx(m0); Vy(m0)];
            Vn0 = [Vx(n0); Vy(n0)];
            
            %angle of collision normal vector and x axis
            CosTheta = dot(ex,normal)/(norm(normal)*norm(ex));
            Thetad = acosd(CosTheta);
            
            %Rotation matrix between n-t and x-y
            Rot = [cosd(Thetad) sign(normal(2))*sind(Thetad);  ...
                -sign(normal(2))*sind(Thetad) cosd(Thetad)];
            % Transforms colliding particles velocities in normal-tangent
            % coordinate
            VVm0 = Rot*Vm0;
            VVn0 = Rot*Vn0;
            % using projected velocities, calculates the after-collicion
            % velocities
            VV_after_m0 = ...
                [((1-e*alpha)*VVm0(1)+alpha*(1+e)*VVn0(1))/(1+alpha);...
                VVm0(2)];
            VV_after_n0 = ...
                [((1+e)*VVm0(1)+(alpha-e)*VVn0(1))/(1+alpha); VVn0(2)];
            % Transforms back the after-collisio velocities to X-Y
            % coordinates
            V1 = Rot'*VV_after_m0;
            V2 = Rot'*VV_after_n0;
            
            % gets after-contact velocity components
            Vx(m0,:) = V1(1);
            Vy(m0,:) = V1(2);
            Vx(n0,:) = V2(1);
            Vy(n0,:) = V2(2);
        end
        %------------------------------------------------------------------
        % updating H matrix which contains the minimum time-to hits between
        % any 2 particles, Note: updating the row and colums
        % associated to the collisded particles by doing calculations,
        % Note: the rest of time-to-hits will be subtracted by previous
        % minimum time-to-hit
        % Note: H(m,n) = inf means the 2 particles m and n never will
        % collide
        H = H - t_hitParticle_min;
        H(H==0) = inf;
        for mm = min(m0,n0):max(m0,n0)
            for nn = setdiff(mm:NumOfParticles ,mm)
                % calculating relative distances and velocities
                dx =  X(mm,:)- X(nn,:);
                dy =  Y(mm,:)- Y(nn,:);
                dVx = Vx(mm) - Vx(nn);
                dVy = Vy(mm) - Vy(nn);
                
                dR = [dx dy];
                dV = [dVx dVy];
                
                if dot(dV,dR) >= 0
                    H(mm,nn) = inf;
                else
                    droot = (dot(dV,dR))^2-dot(dV,dV)*(dot(dR,dR)...
                        -(radii(mm)+radii(nn))^2);
                    if droot <0; H(mm,nn) = inf;
                    else
                        % time-to-hit is calculated here
                        H(mm,nn) = (-dot(dV,dR)-sqrt(droot))/dot(dV,dV);
                    end
                    if dot(dV,dV) <= eps; H(mm,nn) = inf;end
                end
                % colliding particle positions and contact-time are updates
                X0(m0,:) = X(m0,:);
                Y0(m0,:) = Y(m0,:);
                X0(n0,:) = X(n0,:);
                Y0(n0,:) = Y(n0,:);
                t_cont(m0) = t;
                t_cont(n0) = t;
            end
        end
        for nn = min(m0,n0):max(m0,n0)
            for mm = setdiff(1:nn,nn)
                dx =  X(mm,:)- X(nn,:);
                dy =  Y(mm,:)- Y(nn,:);
                dVx = Vx(mm) - Vx(nn);
                dVy = Vy(mm) - Vy(nn);
                
                dR = [dx dy];
                dV = [dVx dVy];
                
                if dot(dV,dR) >=0
                    H(mm,nn) = inf;
                else
                    droot = (dot(dV,dR))^2-dot(dV,dV)*(dot(dR,dR)...
                        -(radii(mm)+radii(nn))^2);
                    if droot <0; H(mm,nn) = inf;
                    else
                        % time-to-hit is calculated here
                        H(mm,nn) = (-dot(dV,dR)-sqrt(droot))/dot(dV,dV);
                    end
                    if dot(dV,dV) <= eps; H(mm,nn) = inf;end
                end
            end
        end
        
        t_hitParticle_min = min(min(H));
        [m0,n0] = find(H == t_hitParticle_min);
        if m0 >1; NumOfSimCont = length(m0);end
        t_0 = t_0 + t_hitParticle_min;
    end
    %% ====================================================================
    pause(0*dt(k))                          % pause the motion for dt
    k=k+1;
      dt(k) = min(t_hitParticle_min, ...
          t_hitParticle_min/floor(t_hitParticle_min/dt(1)));  
%     dt(k) = t_hitParticle_min/2^ndivision;  % updates the time step
%     if dt(k) == inf
%         dt(k) = 0.5*min(radii)/(max(max(abs(Vx),max(abs(Vy)))));
%     end
    % ---------------------------------------------------------------------
    %% Trying to obtain smoother dt variations!
%     D0 = t_hitParticle_min/1.1;
%     if D0 ~= 0 && D0 ~= inf
%         timeStepRatio = dt(k)/dt(k-1);
%         if timeStepRatio <1
%             while dt(k) <= D0
%                 dt(k) = dt(k)*2;
%             end
%         end
%         if timeStepRatio >1
%             while dt(k) >= D0
%                 dt(k) = dt(k)/2;
%             end
%         end
%     end
    %   ===================================================================
    t = t+dt(k);                          % updates time 
end
%delets waitbar
% delete(WaitBar); 
close all