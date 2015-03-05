function dy = self_check(t,y)
%% check ode system result
%% check question 4
dy = zeros(2,1);    % a column vector
dy(1) = 2*y(1)-y(2)^2+t*exp(t)-t;
dy(2) = y(1);

% options = odeset('RelTol',1e-8,'AbsTol',[1e-8 1e-8]);
% [T,Y] = ode45(@rigid,[0 10],[1 0],options);
% plot(T,Y(:,1),'-',T,Y(:,2),'.')

%% 
% %% check question 6
% dy = zeros(2,1);    % a column vector
% dy(1) = -4*y(1)+3*y(2)+6*sin(t^2);
% dy(2) = -2.4*y(1)+1.6*y(2)+3.6*cos(t^2);
% 
% % options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);
% % [T,Y] = ode45(@rigid,[0 4],[0 1],options);
% % plot(T,Y(:,1),'-',T,Y(:,2),'.')