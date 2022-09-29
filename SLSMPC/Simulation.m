classdef Simulation
    properties
        sys;
        state_constr;
        input_constr;
        terminal_constr;
    end

    methods
        function obj = Simulation(sys, Xc, Uc, Tc)
            obj.sys = sys;
            obj.state_constr = Xc;
            obj.input_constr = Uc;
            obj.terminal_constr = Tc;
        end


        function plot_trajectory(obj, traj_set)
            
            num_traj = length(traj_set);
            
            figure;
            obj.state_constr.plot('color', 'lightyellow')
            hold on 
            for ii = 1:num_traj
                plot(traj_set{ii}.x_seq(1,:), traj_set{ii}.x_seq(2,:), 's-', 'LineWidth', 1.5);
                hold on
            end
            
            figure;
            for ii = 1:num_traj
                plot(traj_set{ii}.u_seq, 's-', 'LineWidth', 1.5);
                hold on
            end
            xlabel('$time$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18);
            %title([title_str, ' control inputs, horizon = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);
        end

        function traj_set = simulate(obj, is_open_loop, n_step, num_traj, init_x_minus, init_u_minus, x0, mpc)
            traj_set = cell(1, num_traj);
            
            nu = obj.sys.nu;
            nx = obj.sys.nx;
            epsA = obj.sys.epsA;
            epsB = obj.sys.epsB;
            for ii = 1:num_traj
                traj = struct;
                x = x0;
                x_minus = init_x_minus;
                u_minus = init_u_minus;
               
                % need to flip A because the order of A is {A_0, A_1, ..., A_na} but x_minus: [x_{-na}; ...; x{-1}]
                
                x_seq = [x0];
                u_seq = [];
                
                for jj = 1:n_step
                    A = obj.sys.A;
                    B = obj.sys.B;
                    
                    if isa(mpc, 'PolyAffineTimeDelayMPC')
                       vert1 = mpc.poly_uncertainty_set{1};
                       vert2 = mpc.poly_uncertainty_set{2};
                       lam = rand(1);
                       for i = 1:obj.sys.na+1
                           A{i} = A{i} + lam*vert1.Delta_A{i} + (1-lam)*vert2.Delta_A{i};
                       end
                       for i = 1:obj.sys.nb+1
                           B{i} = B{i} + lam*vert1.Delta_B{i} + (1-lam)*vert2.Delta_B{i};
                       end
                    end

                    
                    A = cell2mat(flip(A));
                    B = cell2mat(flip(B));
                    disp(jj);
                    w = rand(nx,1).*2*obj.sys.sigma_w - obj.sys.sigma_w;
                    if is_open_loop
                        if isa(mpc, 'PolyAffineTimeDelayMPC')
                            u = mpc.K((jj - 1)*nu + 1:jj*nu, 1:jj*nx) * (x_seq - mpc.h(1:jj*nx));
                        else
                            u = mpc.K((jj - 1)*nu + 1:jj*nu, 1:jj*nx) * x_seq;
                        end
                    else
                        try
                            mpc.solve(x, x_minus, u_minus);
                            if  isa(mpc, 'PolyAffineTimeDelayMPC')
                                u = mpc.K(1:nu, 1:nx)*(x - mpc.h(1:nx));
                            else 
                                u = mpc.K(1:nu, 1:nx)*x;
                            end
                        catch
                            break;
                        end
                    end
                    u_seq = [u_seq u];
                    x_plus = A*[x_minus; x] + B*[u_minus; u] + w;
                    x_seq = [x_seq; x_plus];
                    x_minus = [x_minus(nx+1:end); x];
                    if mpc.sys.nb > 0
                        u_minus = [u_minus(nu+1:end); u];
                    end
                    x = x_plus;
                    disp(x);
                end
                traj.x_seq = reshape(x_seq, obj.sys.nx,[]);
                traj.u_seq = reshape(u_seq, obj.sys.nu,[]);
                traj_set{ii} = traj;
              
            end
            
        
        end


        function feasible = isfeasible(obj, traj_set)
            eps = 10^-7;
            feasible = true;
            for i = 1:length(traj_set)
                x_seq = traj_set{i}.x_seq;
                u_seq = traj_set{i}.u_seq;
                state_feasible = all(obj.state_constr.A*x_seq <= obj.state_constr.b + eps, 'all');
                input_feasible = all(obj.input_constr.A*u_seq <= obj.input_constr.b + eps, 'all');
                feasible = state_feasible & input_feasible
            end
        end
        
    end

end
