% This script implements the 3D robust MPC example in 
% Robust Model Predictive Control of Time-Delay Systems through System Level Synthesis
% Conference on Decision and Control, 2022 (https://arxiv.org/abs/2209.11841). 

clear;

nx = 3;
nu = 1;
alpha = 1.2957;
A1 = [
  1.0509 0 0; 
  -0.0509 1 0;
  0.0509*alpha -0.4*alpha 1
];

A1_bar = [
  0.0218 0 0;
  -0.0218 0 0;  
  0.0218*alpha 0 0
];

A_null = zeros(nx, nx);
A = {A1, A_null, A_null, A1_bar}

B = [-0.1429; 0; 0];
B_null = zeros(nx, nu);
B = {B};


% construct Dela_A, Delta_B polytopic set
num_vert = 2;
polytopic_uncertainty_set = cell(1, num_vert);
vert_1 = struct;
vert_1.Delta_A = {[0 0 0; 0 0 0;-0.2957*0.0509 -0.2957*(-0.4) 0], zeros(nx), zeros(nx), [0 0 0; 0 0 0;-0.2957*0.0218 0 0]};
vert_1.Delta_B = {[0; 0; 0]};

vert_2 = struct;
vert_2.Delta_A = {[0 0 0; 0 0 0;0.2957*0.0509 0.2957*(-0.4) 0], zeros(nx), zeros(nx), [0 0 0; 0 0 0;0.2957*0.0218 0 0]};
vert_2.Delta_B = {[0; 0; 0]};

polytopic_uncertainty_set{1}= vert_1; polytopic_uncertainty_set{2} = vert_2;

na = 3;
nb = 0;

% epsA, epsB were used for testing with norm-bounded uncertainty. They are
% unrelated if polytopic model uncertainty is considered. 
epsA = 0; epsB = 0;

sigma_w = 0.05;

time_delay_system = TimeDelaySystem(A, B, epsA, epsB, sigma_w);

Q = eye(nx);
Q_T = 100*Q;
R = 0.01;


%% state and input constraints
Uc = Polyhedron('A', [eye(nu); -eye(nu)], 'b', pi*ones(2*nu,1));

b = [2/3*pi; 2*pi; 15; 4*pi; 2*pi; 15];
Xc = Polyhedron('A', [eye(nx); -eye(nx)], 'b', b);

XT = Xc;

gamma = 1.0;
x0 = gamma*[0.5*pi; 0.75*pi; -5];

x_minus = zeros(nx*na,1);
u_minus = zeros(nu*nb,1);

T = 6;

mpc = PolyAffineTimeDelayMPC(time_delay_system, Xc, Uc, XT, T, Q, Q_T, R);
mpc.assign_poly_uncertainty_set(polytopic_uncertainty_set);

%% closed-loop simulation using robust SLS MPC
simulation = Simulation(time_delay_system, Xc, Uc, XT);
num_traj = 1;
is_open_loop = false;
traj_set = simulation.simulate(is_open_loop, 150, num_traj, x_minus, u_minus, x0, mpc);
%simulation.plot_trajectory(traj_set);
simulation.isfeasible(traj_set)
len = size(traj_set{1}.x_seq, 2); 

%% plot simulated closed-loop trajectories
figure; 
plot(1:len, traj_set{1}.x_seq(1,:),'b', 1:len, traj_set{1}.x_seq(2,:), 'g', ...
    1:len, traj_set{1}.x_seq(3,:),'r', 'LineWidth',1.5);
legend('x1','x2', 'x3', 'FontSize', 14, 'Interpreter', 'Latex');
xlabel('time', 'FontSize', 18, 'Interpreter', 'Latex'); 
ylabel('state', 'FontSize', 18, 'Interpreter', 'Latex');
grid on
xlim([0, 140]);

figure;
plot(1:len-1, traj_set{1}.u_seq(1,:),'b', 'LineWidth',1.5);
legend('u', 'AutoUpdate','off', 'FontSize', 14, 'Interpreter', 'Latex');

hold on 
xlabel('time', 'FontSize', 18, 'Interpreter', 'Latex'); 
ylabel('state', 'FontSize', 18, 'Interpreter', 'Latex');
plot(1:len-1, pi*ones(1, len-1), 'r-.', 'LineWidth', 1.0);
plot(1:len-1, -pi*ones(1, len-1), 'r-.', 'LineWidth', 1.0);
xlabel('time', 'FontSize', 18, 'Interpreter', 'Latex'); 
ylabel('control input', 'FontSize', 18, 'Interpreter', 'Latex');
grid on
xlim([0, 140]);


