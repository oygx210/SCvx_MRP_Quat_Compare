%% Padraig Lysandrou 2020 scvx
% An SCVX 6dof guidance instance
close all; clear all; clc;
iter_limit = 20;
K = 50;
% load vehicle parameters
vehicle_params;

% Configure the initial condition
tf_guess = 10;
r_N_0 = [100 100 100].';
v_N_0 = -[1 1 1].';
sigma = [0 0 0].';
omega = [0.01 0.01 0.01].';

% Can use any form of attitude formalism for sigma here.
lander_nd = compute_nd_factors(lander, tf_guess, r_N_0, v_N_0, sigma, omega);
% Terminal state (terminal mass should be left unconstrained)
lander_nd.XT = zeros(size(lander_nd.X0));

% Compute a prior reference state and control input history over K
[x_0, u_0] = initialize_reference_trajectory(lander_nd, K);

% Set the dimensions of things!
lander_nd.m = length(x_0(:,1));
lander_nd.n = length(u_0(:,1));
lander_nd.K = K;

% Configure weights
w_nu 		= 1.e5;
w_delta 	= 1.e-3;
w_delta_sigma = 0.1; 
weights = [w_nu w_delta w_delta_sigma];

% set up beginning of loop
sigma = tf_guess;
converged = false;
error_flag = false;
iter_counter = 0;
solve_time = 0;
x = lander_nd.X0;
u = lander_nd.U0;

% precompute jacobians for use later
lander_dynamics = vehicle_dynamics(lander_nd);

while ~converged && iter_counter < iter_limit && ~error_flag
   iter_counter = iter_counter + 1;
   disp("Iteration " + string(iter_counter));
   
   % Compute the linear system matrices
   A_bar, B_bar, C_bar, Sigma_bar, z_bar = lander_dynamics.discretized_dynamics(x, u, sigma);
    
   % Form and solve the convex sub-problme
    [x_cvx, u_cvx, sigma_cvx, delta_norm_cvx, ...
        sigma_norm_cvx, nu_cvx, error, stats] = scvx_subproblem(lander, A_bar, ...
                                                B_bar, C_bar, Sigma_bar, ...
                                                z_bar, x, u, sigma, ...
                                                weights, K);
   
    % Flag an error, if that happens
    if error
        error_flag = true;
    else
        x = x_cvx.value();
        
    end
   
end



% while not converged and conv_counter < conv_limit and not error_flag:
%     conv_counter +=1;
%     print("Iteration " + str(conv_counter) + "       CURRENT RUN COUNT " + str(j));
%     # this performs a discrete integration and finds the new ABCSz values
%     # at each time for the next optimization intseration
%     A_bar, B_bar, C_bar, Sigma_bar, z_bar = \
%         integrator.discretized_dynamics(x, u, sigma);
% 
%     # Solves the individual convex problem
%     x_cvx, u_cvx, sigma_cvx, delta_norm_cvx, \
%         sigma_norm_cvx, nu_cvx, error, stats = \
%         cvx_problem(lander, A_bar, B_bar, C_bar, Sigma_bar, \
%             z_bar, x, u, sigma, weights, K);
% 
%     if error:
%         print(error);
%         error_flag = 1;
%     else:
%         x = x_cvx.value;
%         u = u_cvx.value;
%         sigma = sigma_cvx.value;
%         delta_norm = delta_norm_cvx.value;
%         sigma_norm = sigma_norm_cvx.value;
%         nu = nu_cvx.value;
% 
%         nu_norm = np.linalg.norm(nu, np.inf);
% 
%         if delta_norm > 1e-3:
%             weights[1] *= 1.5;
%             # weights[1] *= 20;
%         else:
%             weights[1] *= 1.1
% 
%         if delta_norm < 1e-3 and sigma_norm < 1e-3 and nu_norm < 1e-7:
%             converged = True;
% 
%         solve_time += stats.solve_time + stats.solve_time;
% 
%         thing = np.array([nu_norm, delta_norm, sigma_norm]);
%         print("[nu_norm delta_norm sigma_norm ] ==	" + str(thing));
%         print("Has converged:  		" + str(converged));
%         print("Final Time for Iterate:  " + str(sigma));
%         print("Solve time for Iterate:  " + str(stats.solve_time + stats.solve_time));
%         print("")
%         # store the iteration number and the sigma
%         iterate_data = np.append(iterate_data, [[conv_counter, sigma, x[0,-1]]],  axis=0);
% 
% if not error_flag and converged:
%     solve_times[j] = solve_time;
%     x_total[j*K:(j+1)*K,:] = np.transpose(x);
%     u_total[j*K:(j+1)*K,:] = np.transpose(u);
%     scaling_conts[j,:] = lander.things_that_matter;
%     sigma_total[j] = sigma;
%     j += 1;

