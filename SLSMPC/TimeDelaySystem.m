classdef TimeDelaySystem
   properties
      na; nb; % length of delay
      A; % A is a cell of length na + 1. Ｔhe order is A = {A(0), A(n_a), ..., A(1)} and same for B
      B; % B is a cell of length nb + 1
      epsA = 0; epsB = 0; % magnitudes of model uncertainy
      nx; nu; % dim of state space and input space
      sigma_w = 0; % the norm of additive disturbances w(k)
   end
   
   methods
      function obj = TimeDelaySystem(A, B, epsA, epsB, sigma_w)
        obj.A = A;
        obj.B = B;
        obj.epsA = epsA;
        obj.epsB = epsB;
        obj.sigma_w = sigma_w;
        obj.na = size(obj.A, 2) - 1;
        obj.nb = size(obj.B, 2) - 1;
        obj.nx = size(obj.A{1}, 2);
        obj.nu = size(obj.B{1}, 2);
      end

       
      function [ZA_block, ZB_block, ZA_minus_block, ZB_minus_block] = dynamics_matrix_partition(obj, T)
         % x = ZA x + ZB u + ZA_minus x_minus + ZB_minus u_minus + w
         % This function computes the matrices ZA, ZB, ZA_minus, and ZB_minus
  
         nx = obj.nx; 
         nu = obj.nu; 
         na = obj.na;
         nb = obj.nb;
           
         % construct ZA and A_minus
         A_bar = kron([eye(T) zeros(T, na + 1)], obj.A{1 + na});
         for k = 1:na
            base = [zeros(T, k) eye(T) zeros(T, na + 1 - k)];
            A_bar = A_bar + kron(base, obj.A{na + 1 -k});
         end
            
         ZA = A_bar(:, na*nx+1:end);
         ZA_minus = A_bar(:,1:na*nx);
           
         % construct ZB and ZB_minus
         B_bar = kron([eye(T) zeros(T, nb + 1)], obj.B{1+nb});
         for k = 1:nb
            base = [zeros(T, k) eye(T) zeros(T, nb + 1 - k)];
            B_bar = B_bar + kron(base, obj.B{nb + 1 -k});
         end
            
         ZB = B_bar(:, nb*nu+1:end);
         ZB_minus = B_bar(:,1:nb*nu);
           
         ZA_block = [zeros(nx, size(ZA, 2)); ZA];
         ZA_minus_block = [zeros(nx, size(ZA_minus, 2)); ZA_minus];
         ZB_block = [zeros(nx, size(ZB, 2)); ZB];
         ZB_minus_block = [zeros(nx, size(ZB_minus, 2)); ZB_minus];
           
      end

      %% obtain the Delta_A, Delta_B related matrices
      function [ZDeltaA_block, ZDeltaB_block, ZDeltaA_minus_block, ZDeltaB_minus_block] = uncertainty_matrix_partition(obj, vert, T)
         % vert.Delta_A = {Delta_A_0, ..., Delta_A_na},  vert.Delta_B = {Delta_B_0, ..., Delta_B_nb}
         % This function computes the matrices ZDeltaA, ZDeltaB, ZDeltaA_minus, and ZDeltaB_minus
  
         nx = obj.nx; 
         nu = obj.nu; 
         na = obj.na;
         nb = obj.nb;
         
         Delta_A = vert.Delta_A; Delta_B = vert.Delta_B;
         % construct ZA and A_minus
         A_bar = kron([eye(T) zeros(T, na + 1)], Delta_A{1 + na});
         for k = 1:na
            base = [zeros(T, k) eye(T) zeros(T, na + 1 - k)];
            A_bar = A_bar + kron(base, Delta_A{na + 1 -k});
         end
            
         ZA = A_bar(:, na*nx+1:end);
         ZA_minus = A_bar(:,1:na*nx);
           
         % construct ZB and ZB_minus
         B_bar = kron([eye(T) zeros(T, nb + 1)], Delta_B{1+nb});
         for k = 1:nb
            base = [zeros(T, k) eye(T) zeros(T, nb + 1 - k)];
            B_bar = B_bar + kron(base, Delta_B{nb + 1 -k});
         end
            
         ZB = B_bar(:, nb*nu+1:end);
         ZB_minus = B_bar(:,1:nb*nu);
           
         ZDeltaA_block = [zeros(nx, size(ZA, 2)); ZA];
         ZDeltaA_minus_block = [zeros(nx, size(ZA_minus, 2)); ZA_minus];
         ZDeltaB_block = [zeros(nx, size(ZB, 2)); ZB];
         ZDeltaB_minus_block = [zeros(nx, size(ZB_minus, 2)); ZB_minus];
           
      end
      
      %% compute the d vector which encodes past information
      function [d] = additive_d(obj, T, x_minus, u_minus)
         [ZA_block, ZB_block, ZA_minus_block, ZB_minus_block] = obj.dynamics_matrix_partition(T);
         d = ZA_minus_block*x_minus + ZB_minus_block*u_minus;
      end
     

      %% stack norm of model uncertainty into a matrix
      function [Delta_A_norm, Delta_A_minus_norm, Delta_B_norm, Delta_B_minus_norm] = model_uncertainty_norm(obj, T)
         epsA = obj.epsA(end);
         Delta_A_norm = [epsA*eye(T) zeros(T, obj.na )];
         for k = 1:obj.na
            if size(obj.epsA, 1) ~= 1
                epsA = obj.epsA(end-k);
            end
            Delta_A_norm = Delta_A_norm + [zeros(T, k) epsA*eye(T) zeros(T, obj.na - k)];
            
         end
         
         Delta_A_minus_norm = Delta_A_norm(:,1:obj.na);
         Delta_A_norm = Delta_A_norm(:, obj.na+1:end);
         
         epsB = obj.epsB(end);
         Delta_B_norm = [epsB*eye(T) zeros(T, obj.nb)];
         
         for k = 1:obj.nb
            if size(obj.epsB, 1) ~= 1
                epsB = obj.epsB(end-k);
            end
            Delta_B_norm = Delta_B_norm + [zeros(T, k) epsB*eye(T) zeros(T, obj.nb - k)];
         end
         
         Delta_B_minus_norm = Delta_B_norm(:,1:obj.nb);
         Delta_B_norm = Delta_B_norm(:, obj.nb+1:end);

      end

      function [d_delta_upper_bd] = compute_d_delta_norm(obj, T, x_minus, u_minus,  Delta_A_minus_norm, Delta_B_minus_norm)
         x_minus_norm = [];
         for i = 1:obj.na
            x_minus_norm = [x_minus_norm; norm(x_minus(obj.nx*(i-1)+1:obj.nx*i), inf)];
         end
            
         u_minus_norm = [];
         for i = 1:obj.nb
            u_minus_norm = [u_minus_norm; norm(u_minus(obj.nu*(i-1)+1:obj.nu*i), inf)];
         end

         d_delta_upper_bd = Delta_A_minus_norm * x_minus_norm + Delta_B_minus_norm * u_minus_norm;


      end



   end

end