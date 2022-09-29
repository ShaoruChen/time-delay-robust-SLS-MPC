% This script implements the solver time evaluation in 
% Robust Model Predictive Control of Time-Delay Systems through System Level Synthesis
% Conference on Decision and Control, 2022 (https://arxiv.org/abs/2209.11841). 

clear;
na_list = [8 16 24 32 40];
num_list = [20 20 20 20 20];

%% solver time comparison
results = cell(1, 5);
for ii = 1:length(na_list)
    ii
    
    na = na_list(ii);
    nb = na/2;
    T = na + 5;
    num = num_list(ii);
    
    [sol_cell] = solver_time(na, nb, T, num);
        
    results{ii} = sol_cell;
    save solver_time_comparison_temp_T results
end

save solver_time_comparison_T

%% simualtions that do not consider time delays
T_list = [8 16 24 32 40] + 5;
num_list = [20 20 20 20 20];

results_non_delay = cell(1, 5);
for ii = 1:length(na_list)
    ii
    
    na = 0;
    nb = 0;
    T = T_list(ii);
    num = num_list(ii);
    
    [sol_cell] = solver_time(na, nb, T, num);
        
    results_non_delay{ii} = sol_cell;
    save solver_time_comparison_non_delay_temp_T results_non_delay
end
save solver_time_comparison_non_delay_T 

% process data
mean_solver_time_list = zeros(1,5);
for ii = 1:5
    status = zeros(1, 20);
    solver_time = zeros(1, 20);
    for jj = 1:20
        status(jj) = results{ii}{jj}.solution.problem;
        solver_time(jj) = results{ii}{jj}.solver_time;
    end
    valid_solver_time = solver_time(status == 0);
    mean_solver_time = mean(valid_solver_time);
    mean_solver_time_list(ii) = mean_solver_time;
end

mean_solver_time_list_non_delay = zeros(1,5);
for ii = 1:5
    status = zeros(1, 20);
    solver_time = zeros(1, 20);
    for jj = 1:20
        status(jj) = results_non_delay{ii}{jj}.solution.problem;
        solver_time(jj) = results_non_delay{ii}{jj}.solver_time;
    end
    valid_solver_time = solver_time(status == 0);
    mean_solver_time = mean(valid_solver_time);
    mean_solver_time_list_non_delay(ii) = mean_solver_time;
end


