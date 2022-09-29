
classdef PolyAffineTimeDelayMPC < ModelPredictiveControl
    % robust SLS MPC method for time delay systems
    properties (SetAccess = public)
        poly_uncertainty_set =[]; % polytopic uncertainty set vertices 
        num_vert;
        h; 
        d;
    end
    
   
    methods (Access = public)
        function assign_poly_uncertainty_set(obj, poly_uncertainty_set)
            obj.poly_uncertainty_set = poly_uncertainty_set;
            obj.num_vert = length(poly_uncertainty_set);
        end
        
        function [sol] = solve(obj, x0, x_minus, u_minus)
            
            T = obj.T;
            nx = obj.sys.nx;
            nu = obj.sys.nu;
            sigma_w = obj.sys.sigma_w;
            
            Phi_x = sdpvar( (T + 1) * nx, (T + 1) * nx, 'full');
            Phi_u = sdpvar( (T + 1) * nu, (T + 1) * nx, 'full');
            
            % construct the sigma matrix
            sigma_seq = sdpvar(1, T*nx);
            
            % Manually set the off-diagonal entries as zero. The full parameterization 
            % of the filter Sigma_mat can be applied too (see the commented
            % lines).
            Sigma_mat = eye(nx);
            Sigma_mat = blkdiag(Sigma_mat, diag(sigma_seq));
%             Sigma_mat = sdpvar((T+1)*nx, (T+1)*nx, 'full');

            % compute h
            [ZA_block, ZB_block, ZA_minus_block, ZB_minus_block] = obj.sys.dynamics_matrix_partition(T);
            Id = eye((T + 1)*nx);
            d = ZA_minus_block*x_minus + ZB_minus_block*u_minus;
            obj.d = d;
            h = inv(Id - ZA_block)*d;
            

            % The cost, no penalty on the initial state and u_T
            % nominal_x is [x_1, x_2, ... ,x_T]
            % nominal_u is [u_0, u_1, ... ,u_{T-1}]
            Q_bar = kron(eye(T+1), obj.Q);
            Q_bar(T*nx+1: (T+1) * nx, T*nx+1:(T+1) * nx) = obj.Q_T;
            for ii = 1:obj.sys.na
                Q_bar((T-ii)*nx+1: (T-ii+1) * nx, (T-ii)*nx+1:(T-ii+1) * nx) = obj.Q_T;
            end
            
            R_bar = kron(eye(T+1), obj.R);
            nominal_x = Phi_x(1:nx*(T+1), 1:nx) * x0 + h;
            nominal_u = Phi_u(1:nu*(T+1), 1:nx) * x0;
            cost = (nominal_x')*Q_bar*nominal_x + (nominal_u')*R_bar*nominal_u;


            % add constraints to the problem
            constr = [];
            
%           % structure constraint on Sigma_mat if the full parameterization
%           % is used.
% 			for k = 1 : T
% 				constr = [constr, Sigma_mat( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (T + 1 - k)*nx)];
%             end
% 
%             constr = [constr, Sigma_mat(1:nx, 1:nx) == eye(nx)];
%             for k = 1:T
%                 constr = [constr, Sigma_mat(k*nx+1:(k+1)*nx, k*nx+1:(k+1)*nx) == diag(sigma_seq((k-1)*nx+1:k*nx))];
%             end
            
            % Phi_x and Phi_u are block lower triangular matrix
            for k = 1 : T
                constr = [constr, Phi_x( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (T + 1 - k)*nx)];
            end

            for k = 1 : T
                constr = [constr, Phi_u( (k - 1)*nu + 1: k*nu, k*nx + 1: end) == zeros(nu, (T + 1 - k)*nx)];
            end
            
            % add affine constraint
            constr = [constr, (Id - ZA_block)*Phi_x - ZB_block*Phi_u == Sigma_mat];
            
            % state constraints
            if ~isempty(obj.state_constr)
                Fx = obj.state_constr.A; 
                bx = obj.state_constr.b;
                nFx = size(Fx, 1); 
                nbx = length(bx);
            else
                warning('must have state constraints');
            end
            
            for ii = 1:T
                for jj = 1: nFx
                    f = Fx(jj,:); b = bx(jj);
                    LHS = f*Phi_x((ii-1)*nx+1:ii*nx,1:nx)*x0;
                    for kk = 1:ii-1
                        LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx),1);
                    end
                    LHS = LHS + f*h((ii-1)*nx+1:ii*nx);
                    constr = [constr, LHS <= b];  
               end
            end

            % terminal_constr
            if ~isempty(obj.terminal_constr)
                Ft = obj.terminal_constr.A;
                bt = obj.terminal_constr.b;
                nFt = size(Ft, 1);
            else
                Ft = Fx;
                bt = bx;
                nFt = nFx;
            end
            
            % add terminal constraints only on the last state
%             for jj = 1: nFt
%                 f = Ft(jj,:); b = bt(jj);
%                 LHS = f*Phi_x( T*nx+1: (T+1)*nx,1:nx)*x0;
%                 for kk = 1:T
%                     LHS = LHS + norm(f*Phi_x(T*nx+1: end, kk*nx+1: (kk+1)*nx),1);
%                 end
%                 LHS = LHS + f*h(T*nx+1:(T+1)*nx);
%                 constr = [constr, LHS <= b]; 
%             end

            % add terminal constraints on the last n_a states
            for ii = T - obj.sys.na+1:T+1
                for jj = 1: nFt
                    f = Ft(jj,:); b = bt(jj);
                    LHS = f*Phi_x((ii-1)*nx+1:ii*nx,1:nx)*x0;
                    for kk = 1:ii-1
                        LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx),1);
                    end
                    LHS = LHS + f*h((ii-1)*nx+1:ii*nx);
                    constr = [constr, LHS <= b];  
               end
            end
            
            % add input constraint
            if ~isempty(obj.input_constr)
                Fu = obj.input_constr.A; 
                bu = obj.input_constr.b;
                nFu = size(Fu, 1); 
                nbu = length(bu);
            else
                warning('must have input constraints');
            end
            
            for ii = 1:T+1
                for jj = 1: nFu
                    f = Fu(jj,:); b = bu(jj);
                    LHS = f*Phi_u((ii-1)*nu+1:ii*nu,1:nx)*x0;
                    for kk = 1:ii-1
                        LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),1);
                    end
                    constr = [constr, LHS <= b];  
               end
            end

            % uncertainty set over-approximation constraints
            num_vert = obj.num_vert;
            eye_mat = eye(nx);
            
%             % if full parameterization of Sigma_mat is applied:  
%             Sigma_diag = diag([ones(1,nx) sigma_seq]);
%             Sigma_subdiag = Sigma_mat - Sigma_diag;

            % We have the following equation holds:
            % Sigma_diag tilde(w) = (ZDeltaA*Phi_x + ZDeltaB*Phi_u -
            % Sigma_subdiag)*tilde(w) + ZDeltaA*h + ZDeltaA_minus x_minus
            % + ZDeltaB_minus u_minus + w
            for ii = 1:num_vert
               vert = obj.poly_uncertainty_set{ii};
               [ZDeltaA_block, ZDeltaB_block, ZDeltaA_minus_block, ZDeltaB_minus_block] = ...
                                                        obj.sys.uncertainty_matrix_partition(vert, T);

               block_A = ZDeltaA_block*Phi_x + ZDeltaB_block*Phi_u;
                               
%              % if full parameterization of Sigma_mat is applied:                                      
%              block_A = ZDeltaA_block*Phi_x + ZDeltaB_block*Phi_u - Sigma_subdiag;

               block_known = block_A(:, 1:nx)*x0 + ZDeltaA_block*h + ZDeltaA_minus_block*x_minus + ...
                                                                     ZDeltaB_minus_block*u_minus;
               block_unknown = block_A(:,nx+1:end);
               for jj = 1:T
                   for kk = 1:nx
                       e_vec = eye_mat(kk, :);
                       temp = norm(e_vec*block_known(jj*nx+1:(jj+1)*nx), inf);
                       for tt = 0:jj-1
                           temp = temp + norm(e_vec*block_unknown(jj*nx+1:(jj+1)*nx, tt*nx+1:(tt+1)*nx),1);
                       end
                       temp = temp + sigma_w;
                       constr = [constr, temp <= sigma_seq((jj-1)*nx+kk)];
                   end
               end
            end
            
            
            ops = sdpsettings('solver', 'mosek', 'verbose', 2);
            solution = optimize(constr, cost, ops);

            if solution.problem ~= 0
                warning(['Numerical error detected. Error code ', num2str(solution.problem), '\n']);
            end

            sol = struct;
            Phi_x_val = value(Phi_x); 
            Phi_u_val = value(Phi_u);
            Phi_u_val(find(isnan(Phi_u_val)==1)) = 0;
            fval = value(cost);
            sigma_seq_val = value(sigma_seq);
            
            sol.Phi_x = Phi_x_val; sol.Phi_u = Phi_u_val;
            sol.fval = fval; 
            sol.sigma_seq = sigma_seq_val;
            sol.status = solution.problem;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;
            sol.d = d;
            obj.h = h;
            obj.K = sol.Phi_u*inv(sol.Phi_x);
            yalmip('clear');

        end
    end
end