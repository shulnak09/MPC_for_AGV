close all
clear 
clc



% Define the state variables
syms X Y psi vx vy r delta F_x d_delta

% Generate Desired trajectory using Waypoint Generator:
wp = 1*[0 0 0; 198.599 -89.868 0;  401.9228 -181.6857 0; 573.4097 -266.6172  0; 698.5358 -353.4334 0];


%% Generating a Astroid Trajectory
a = 20;
% t = linspace(-1*pi,1*pi,100);
% x = round(a*sin(5*t).*cos(t),2);
% y = round(a*sin(5*t).*sin(t),2);
% x = a * cos(t).^3;
% y = a * sin(t).^3;
% % 
% wp = [x',y',zeros(100,1)];


% time_instance = [0.1  0.2  0.8];
% c = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980]};
% fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);
% for i  = 1:length(time_instance)
toa = 4*(0:size(wp,1)-1).';
Fs = 20;




[pos,ori,vel_x,vel_y,yaw_rate,steer] = trajectory_generator(Fs,toa,wp);
x_ref = [pos,ori,vel_x,vel_y,yaw_rate,steer]';

% figure('Renderer', 'painters', 'Position', [10 10 600 400])
% 
% plot(x_ref(1,:),x_ref(2,:),'b','LineWidth',2)
% ylabel('Y')
% xlabel('X')
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% 
% f = gcf;
% exportgraphics(f,'Trajectory_8_fiure.png','Resolution',300)
% exportgraphics(f,'Trajectory_8_figure.pdf','ContentType','vector',...
%                'BackgroundColor','none')

x  =[X; Y; psi; vx; vy; r; delta];
u = [F_x
    d_delta];
% length_values = [0.18 * 0.5, 0.18, 0.18 *2, 0.18 *3];
% Define the constants:
% for i  = 1:length(length_values)
m = 1500;      % mass of the RC car
I_z = 1500;  % Moment of Inertia
l_r = 1.35;    % length of the rear axle from centre
l_f = 1.62 ;    % length of front axle from centre
g = 9.81;      % acceleration due to gravity
B = 7.4;       % Pajecka Tire constants
C = 1.25;
mu = 0.234;
F_z = m * g/2; % Normal force on each tire
D = mu * F_z; 

% Define the slip angles:
alpha_R = atan((x(5) - l_r *x(6))/x(4));
alpha_F = atan((x(5) + l_f *x(6))/x(4)) - x(7);

F_Ry = D * sin(C*atan(B*alpha_R));
F_Fy = D * sin(C*atan(B*alpha_F));

% % % Define the continuous state space:
f_c = [x(4)*cos(x(3)) - x(5)*sin(x(3))
       x(4)*sin(x(3)) + x(5)*cos(x(3))
       x(6)
       u(1)/m
       (x(4)*u(2) + x(7)*u(1)/m)*(l_r/(l_r+l_f))
       (x(4)*u(2) + x(7)*u(1)/m)*(1/(l_r+l_f))
       u(2)
        ];




% create matlab function for the symbolic expressios:
x_func = matlabFunction(x, 'File', 'x', 'vars', x);
u_func = matlabFunction(u, 'File', 'u', 'vars', u);
f_c_func = matlabFunction(f_c, 'File', 'f_c', 'vars', {x,u});


% Disceretizn the state space based on Euler approximation:
% x(t+1) = x(t) + Ts * f_c 
Ts = 1/Fs;
f_d = x + Ts * f_c;
f_d_func = matlabFunction(f_d, 'File', 'f_d', 'vars', {x,u});

% Compute the jacobian of the state space representation:
% Discrete state matrices:
Ad = jacobian(f_d, x);
Bd = jacobian(f_d, u);

Ad_func = matlabFunction(Ad, 'File', 'Ad', 'vars', {x,u});
Bd_func = matlabFunction(Bd, 'File', 'Bd', 'vars', {x,u});



% compute the equilibrium point for the eqn : x_dot = f_d(x,u)
sol = vpasolve(f_d ==0,3);
x_eq = double([sol.X,sol.Y,sol.psi,sol.vx,sol.vy,sol.r,sol.delta]);
u_eq = double([sol.F_x,sol.d_delta]);
x_eq = x_eq';
u_eq = u_eq';


% Compute the discrete state A_d and B_d, d:
A_d = Ad_func(x_ref(:,1),zeros(2,1));
B_d = Bd_func(x_ref(:,1),zeros(2,1));




% Apply MPC to the original Non-linear system:
N = 10;                 % control Horizon
P = 1e-3* diag([1e6;1e6;1e3;2e5;2e5;1e3;1e3]);           % Terminal State Penalty 'P'
Q = 1e-3* diag([1e6;1e6;1e3;2e5;2e5;1e3;1e3]);        % State Penalty
R = diag([1e1;1e1]);        % Control Penalty
N_sample = size(pos,1);   % Total Samples

% Initial State
% x0 = [0.05; -0.05; 0.1; 0.1; 0.4; 0.25; 0.3];
x0 = x_ref(:,1);
d = x0 - A_d * x0 - B_d * u_eq;

% A_d_d = A_d * Ts + eye(7);
% B_d_d = B_d * Ts;
% d_d = d * Ts;

[n,~] = size(A_d);
[~,m] = size(B_d);

% Define the Upper and lower bounds forming the polyhedra:
x_lb = [-1000; -1000; -10*pi; -100; -20; -2*pi; -0.5*pi];
x_ub = [ 1000;  1000; +40*pi; +100; +20; +2*pi; +0.5*pi];

u_lb = [-20000; -1*pi];
u_ub = [+20000; +1*pi];

xf_lb = [-1000; -1000; -10*pi;  -100; -20; -2*pi; -0.5*pi];
xf_ub = [ 1000;  1000; +40*pi;  +100; +20; +2*pi; +0.5*pi];


% Evolution of state space: 
x_state_evol = zeros(n,N_sample);
x_state_evol(:,1) = x0;
% f_state_evol(:,1) = f_d;

% Evolution of control input:
u_state_evol =zeros(m, N_sample);

% Time
Time = (0:1:N_sample-1);
t = 0; % Simulation time



for k = 1:N_sample-1

    x_curr = x_state_evol(:,k);  % Measure the real state


%   
%     if  (k< N_sample - N)
%         x_des  = x_ref(:,k:k+N-1);
%     else
%         x_des = x_ref(:,N_sample-N+1:N_sample);
%     end
    
    x_des = x_ref(:,k);
    % Prepare the parameters like equality and inequality constraints for
    [H, f, Aeq, beq, lb, ub] = Prepare_QP_solver(A_d, B_d, d, x_lb, x_ub, xf_lb, xf_ub, u_lb, u_ub, P, Q, R, N, x_curr, x_des);
    
   
    % Compute the optimal state solution using Quadratic Programming:
    options = optimoptions('quadprog','Display', 'off');

    % Solve for optimal control using QP solver
    [Z_star, f_star, exitflag, output , lambda] = quadprog(H, f, [], [], Aeq, beq, lb, ub, x_curr, options);

%          disp(Z_star)
    
    
    u_MPC = Z_star(n*N+1:n*N+m);
    u_state_evol(:,k) = u_MPC;
    
%     
% 
%     disp(f_c_init)

    % Evolve the Non-Linear State:
    [t_ode,state] = ode45(@(t,x,u) non_linear_dynamics(t,x,u_MPC,m,I_z,l_r,l_f,g,B,C,mu),[t,t+Ts],x_curr);
    t = t + Ts;
    x_state_evol(:, k+1) = state(end,:)';

    % Compute the discrete state A_d and B_d, d:
    A_d = Ad_func(x_state_evol(:, k+1), u_state_evol(:,k));
    B_d = Bd_func(x_state_evol(:, k+1), u_state_evol(:,k));

    f_d_state = f_d_func(x_curr,u_MPC);
    d = f_d_state - A_d * x_curr - B_d *  u_state_evol(:,k);

end

samples = size(x_state_evol,2);



figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(4,2,1)
plot(Time(1:samples)*Ts,x_state_evol(1,1:samples),'LineWidth',2)
hold on
plot(Time(1:samples)*Ts,x_ref(1,1:samples),'LineWidth',2)
legend('MPC','des')
ylabel('x')

subplot(4,2,2)
plot(Time(1:samples)*Ts,x_state_evol(2,1:samples),'LineWidth',2)
hold on
plot(Time(1:samples)*Ts,x_ref(2,1:samples),'LineWidth',2)
ylabel('y')

subplot(4,2,3)
plot(Time(1:samples)*Ts,x_state_evol(3,1:samples),'LineWidth',2)
hold on
plot(Time(1:samples)*Ts,x_ref(3,1:samples),'LineWidth',2)
ylabel('\psi')

subplot(4,2,4)
plot(Time(1:samples)*Ts,x_state_evol(4,1:samples),'LineWidth',2)
hold on
plot(Time(1:samples)*Ts,x_ref(4,1:samples),'LineWidth',2)
ylabel('v_{x}')

subplot(4,2,5)
plot(Time(1:samples)*Ts,x_state_evol(5,1:samples),'LineWidth',2)
hold on
plot(Time(1:samples)*Ts,x_ref(5,1:samples),'LineWidth',2)
ylabel('v_{y}')


subplot(4,2,6)
plot(Time(1:samples)*Ts,x_state_evol(6,1:samples),'LineWidth',2)
hold on
plot(Time(1:samples)*Ts,x_ref(6,1:samples),'LineWidth',2)
ylabel('r')

subplot(4,2,7)
plot(Time(1:samples)*Ts,x_state_evol(7,1:samples),'LineWidth',2)
hold on
plot(Time(1:samples)*Ts,x_ref(7,1:samples),'LineWidth',2)
ylabel('\delta')
xlabel('Time(s)')
hold off
set(findall(gcf,'-property','FontSize'),'FontSize',14)
f = gcf;
exportgraphics(f,'Figure8.png','Resolution',300)
exportgraphics(f,'Figure8.pdf','ContentType','vector',...
               'BackgroundColor','none')


figure('Renderer', 'painters', 'Position', [10 10 900 600])
f = cell(size(Time,2),1);
for i  = 1:size(Time,2)
    plot(x_state_evol(1,:),x_state_evol(2,:), 'LineWidth',2, Color='b')
    plot(x_ref(1,:),x_ref(2,:),'LineWidth',2)
    xlim([-0,1000])
    ylim([-500,0])
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    hold on
    L = 30; % width of car
    W = 18;  % length of car
    xCenter = x_state_evol(1,i); 
    yCenter = x_state_evol(2,i); 
    
    theta = x_state_evol(3,i);
    cth = cos(theta);
    sth = sin(theta);
    X = [-L/2  L/2  L/2  -L/2];
    Y = [-W/2  -W/2  W/2  W/2];
    Xrot = X.*cth - Y.*sth;
    Yrot = X.*sth + Y.*cth;

    fill(Xrot + xCenter,  Yrot + yCenter, 'r', 'FaceAlpha',  0.3)
%     rectangle('Position',  [xLeft, yBottom, width, height],  'EdgeColor',  'b', 'FaceColor',  'w',  'LineWidth',  2);
    drawnow();
    f{i} = getframe(gcf);
    pause(0.05)
    hold off
end

vidObj = VideoWriter('Figure8_Tracking.avi');
vidObj.Quality = 100;
open(vidObj);
for i = 1:length(f)
    writeVideo(vidObj,f{i})
end
close(vidObj)


% figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(x_state_evol(1,:),x_state_evol(2,:),'LineWidth',2,color = c{i} )
hold on

plot(x_ref(1,:),x_ref(2,:),'r--','LineWidth',2)
hold off
xlim([-10,55])
ylim([-15,15])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
xlabel('X')
ylabel('Y')

f=gcf;

legend('L=0.09 m','L = 0.18 m','L =0.36 m', 'L = 0.54 m', 'des')
exportgraphics(f,'Figure8_Tracking.png','Resolution',300)
exportgraphics(f,'Figure8_Tracking.pdf','ContentType','vector',...
               'BackgroundColor','none')


% 
samples = size(u_state_evol,2);
figure('Renderer', 'painters', 'Position', [10 10 600 400])
subplot(2,1,1)
plot(Time(1:samples)*Ts,u_state_evol(1,1:samples),'LineWidth',2)
ylabel('u(t)')
subplot(2,1,2)
plot(Time(1:samples)*Ts,u_state_evol(2,1:samples),'LineWidth',2)
ylabel('\nabla\delta')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
xlabel('Time(s)')
f = gcf;
exportgraphics(f,'Figure8_control.png','Resolution',300)
exportgraphics(f,'Figure8_control.pdf','ContentType','vector',...
               'BackgroundColor','none')



% 
% figure(3)
% set(gcf, 'PaperPosition', [0 0 4 2]);
% subplot(4,2,1)
% plot(linspace(0,toa(end),size(x_ref,2)),x_ref(1,:))
% legend('x')
% 
% subplot(4,2,2)
% plot(linspace(0,toa(end),size(x_ref,2)),x_ref(2,:))
%  legend('y')
% 
% subplot(4,2,3)
% plot(linspace(0,toa(end),size(x_ref,2)),x_ref(3,:))
% legend('\psi')
% 
% subplot(4,2,4)
% plot(linspace(0,toa(end),size(x_ref,2)),x_ref(4,:))
% legend('v_{x}')
% 
% subplot(4,2,5)
% plot(linspace(0,toa(end),size(x_ref,2)),x_ref(5,:))
% legend('v_{y}')
% 
% subplot(4,2,6)
% plot(linspace(0,toa(end),size(x_ref,2)),x_ref(6,:))
% legend('r')
% 
% subplot(4,2,7)
% plot(linspace(0,toa(end),size(x_ref,2)),x_ref(7,:))
% legend('\delta')
% xlabel('Time')
% 
% subplot(4,2,8)

%  end

