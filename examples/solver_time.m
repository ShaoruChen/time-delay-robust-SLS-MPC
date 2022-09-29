function [sol_cell] = solver_time(na, nb, T, num)
%solver_time: record the solver time of robust SLS MPC on num randomly
%generated discrete-time time-delay systems with delay horizon (na, nb),
%MPC prediction horizon T. 

nx = 2;
nu = 1;

sol_cell = cell(1, num);

for ii = 1:num

A = mat2cell(0.3*randn(nx, nx*(na+1)), [nx], nx*ones(1, na+1));
B = mat2cell(0.3*randn(nx, nu*(nb+1)), [nx], nu*ones(1, nb+1));

% construct polytopic uncertainty set 
eps_A = 0.1;

num_vert = 2;
polytopic_uncertainty_set = cell(1, num_vert);

vert_1 = struct;
Delta_A_1 = repmat([eps_A 0; 0 0], 1, na+1);
Delta_B_1 = repmat(zeros(2,1), 1, nb+1);
vert_1.Delta_A = mat2cell(Delta_A_1, [nx], nx*ones(1,na+1)); vert_1.Delta_B = mat2cell(Delta_B_1, [nx], nu*ones(1,nb+1));

vert_2 = struct;
Delta_A_2 = repmat([-eps_A 0; 0 0], 1, na+1);
Delta_B_2 = repmat(zeros(2,1), 1, nb+1);
vert_2.Delta_A = mat2cell(Delta_A_2, [nx], nx*ones(1,na+1)); vert_2.Delta_B = mat2cell(Delta_B_2, [nx], nu*ones(1,nb+1));
polytopic_uncertainty_set{1} = vert_1; polytopic_uncertainty_set{2} = vert_2;

% construct time delay system
epsA = 0;
epsB = 0;
sigma_w = 0;

time_delay_system = TimeDelaySystem(A, B, epsA, epsB, sigma_w);

Q = eye(nx);
Q_T = Q;
R = 1;

%% state and input constraints
Uc = Polyhedron('A', [eye(nu); -eye(nu)], 'b', 5*ones(2*nu,1));

b = 3*[10; 10; 10 ; 10];
Xc = Polyhedron('A', [eye(nx); -eye(nx)], 'b', b);

XT = Xc;

gamma = 0.5;
x0 = gamma*[5.0; -5.0];
% x0 = [0; 0 ; 0];

x_minus = zeros(nx*na,1);
u_minus = zeros(nu*nb,1);

mpc = PolyAffineTimeDelayMPC(time_delay_system, Xc, Uc, XT, T, Q, Q_T, R);
mpc.assign_poly_uncertainty_set(polytopic_uncertainty_set);
sol = mpc.solve(x0, x_minus, u_minus);

disp(sol.solution.info);
disp(sol.solver_time);

sol_cell{ii} = sol;

end

end

