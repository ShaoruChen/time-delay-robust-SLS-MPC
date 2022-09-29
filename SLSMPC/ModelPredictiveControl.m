classdef ModelPredictiveControl < handle
    
    properties (SetAccess = public)
        sys; % time delay linear system with disturbance
        T;  % horizon used for MPC
        Q; R; Q_T; % quadratic stage cost
        state_constr; % state constraints
        input_constr; % input constraints
        terminal_constr; 
        K; % state feedback controller
    end
    
    methods (Access = public)
        function obj = ModelPredictiveControl(sys, Xc, Uc, Tc, T, Q, Q_T, R)
            obj.sys = sys;
            obj.state_constr = Xc;
            obj.input_constr = Uc;
            obj.terminal_constr = Tc;
            obj.T = T;
            obj.Q = Q; 
            obj.Q_T = Q_T;
            obj.R = R;
        end

    end

    methods (Abstract)
        solve(obj) 
    end

end