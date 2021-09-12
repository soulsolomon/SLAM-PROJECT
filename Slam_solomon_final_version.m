clc
myapp = SPRob;
x = 20;
y = 20;
theta = 0;
speed = 0;
rotationspeed = -0.5;
cov = ones(18,18)*0.00001;
i=0;
myanchorsunit = [-1,-1,-1];
mu = [x,y,theta,myanchorsunit,myanchorsunit,myanchorsunit,myanchorsunit,myanchorsunit]';
dt = 0.1;
myapp.dt=dt;
while(true)
speed = 3.0;
myapp.setspeed(speed);
myapp.setrotationspeed(rotationspeed);
signal = myapp.getanchorssignal();


[mu, cov] = EKFslam(mu,cov,speed,rotationspeed,signal,dt);%Develop this function
% [G] = EKFslam(mu,cov,speed,rotationspeed,dt)

myapp.updateyourAGV(mu(1,1),mu(2,1),mu(3,1));
myanchors=[mu(4,1),mu(5,1);...
mu(7,1),mu(8,1);...
mu(10,1),mu(11,1);...
mu(13,1),mu(14,1);...
mu(16,1),mu(17,1)];
myapp.updateyouranchors(myanchors);

myapp.update;
pause(dt);
end

function [mu_new, cov_new] = EKFslam(mu,cov,speed,rotationspeed,signal,dt);
%-----------------------SETUP-----------------------------------------
%2d plane, translation velocity and rotational velocity which controls 
%the robots movment, speed of rotation, and speed of
%translation, point landmarks (anchors signals), 
%   EKF steps:
%       Initialization:
%       Prediction:
%       Correction:      
%---------------------------------------------------------------------

vt      = speed;        %translational velocity
wt      = rotationspeed;%rotational velocity
theta   = mu(3);        %theta
x       = mu(1);        %pose at x
y       = mu(2);        %pose at y
R       = 1e-19*18;

u       = [signal(:,1) signal(:,2)]; % range-bearing observation

% %eye(3,18) just for mapping to 2n+3 dimension where n is...
% %number of anchors
% update State vector or mu

%---------------Prediction steps--------------------------------------------%%%
% update covariance. We need Jacobian Gtx
% To keep mu uptodate, h() and g() is needed

mu_new  =  mu + eye(3,18)' * [-(vt/wt)*sin(theta) + vt/wt * sin(theta + wt*dt);...
           vt/wt*cos(theta) - vt/wt * cos(theta + wt*dt);...
           wt*dt]

Gxt     = [1 0 -(vt/wt)*cos(theta) + vt/wt * cos(theta + wt*dt);
           0 1 -(vt/wt)*sin(theta) + vt/wt * sin(theta + wt*dt);
           0 0                    1                         ]; %Jacobian
Gt      = eye(18,18) + eye(3,18)' * Gxt * eye(3,18);

%update covariance
cov_new     = Gt * cov * Gt' +  R;
%----------------------------------------------------------------------------%%
%----------------Correction Step and update----------------------------------%%
j = 4;
for i=1:5
   
    landmark_x(i) =  mu_new(1,1) + signal(i,1) * cos( signal(i,2) + mu_new(3,1));
    landmark_y(i) =  mu_new(2,1) + signal(i,1) * sin( signal(i,2) + mu_new(3,1));
    landmark_z(i) =  signal(i,3);
        
    mu_new(j,1)   = landmark_x(i);
    mu_new(j+1,1) = landmark_y(i);
    mu_new(j+2,1) = landmark_z(i);
    j = j + 3;
          
end

plot(mu_new(1),mu_new(2),'b--*')
hold on
title('Robot pose Trajectory');
legend('R_pose predicted value');
xlabel('X');
ylabel('Y');

%--------------------------------------------------------------------------------%%
%--------------------------------END---------------------------------------------%%

end
