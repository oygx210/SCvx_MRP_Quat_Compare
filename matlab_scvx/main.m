%% Padraig Lysandrou 2020 scvx
% An SCVX 6dof guidance instance
close all; clear all; clc;
iter_limit = 20;
K = 30;
% load vehicle parameters
vehicle_params;

% Configure the initial condition
tf_guess = 12.14;
r_N_0 = [0 0 800].';
v_N_0 = -[10 10 80].';
sigmaBN = [0 0 0].';
omega = [0 0 0].';

% Can use any form of attitude formalism for sigma here.
lander_nd = compute_nd_factors(lander, tf_guess, r_N_0, v_N_0, sigmaBN, omega);

% Terminal state (terminal mass should be left unconstrained)
lander_nd.XT = [lander_nd.m_dry; zeros(12,1)];

% Compute a prior reference state and control input history over K
[x_0, u_0] = initialize_reference_trajectory(lander_nd, K);

% Set the dimensions of things!
lander_nd.m = length(x_0(:,1));
lander_nd.n = length(u_0(:,1));
lander_nd.K = K;

% Configure weights
weights.w_nu	= 1.e7;
weights.w_dxu 	= 1.e-3;
weights.w_ds    = 0.1; 
weights.w_s     = 1;

% set up beginning of loop
eta = tf_guess;
converged = false;
error_flag = false;
iter_counter = 0;
solve_time = 0;
x = x_0;
u = u_0;

% precompute jacobians for use later
lander_dynamics = vehicle_dynamics(lander_nd);

while ~converged && iter_counter < iter_limit
   iter_counter = iter_counter + 1;
   disp("Iteration " + string(iter_counter));
   
   % Compute the linear system matrices
   tic
   output_matrices = lander_dynamics.discretized_dynamics(x, u, eta);
    
   % Form and solve the convex sub-problem. Send in the last trajectory and
   % input history for the trust regions.
   o_cvx = scvx_subproblem(lander_nd, output_matrices, x, u, eta, weights);
   disp('dxu_norm   = '+string(o_cvx.delta_norm));
   disp('sigma_norm = '+string(o_cvx.sigma_norm));
   disp('nu_norm    = '+string(o_cvx.nu_norm));
   toc
   
   x   = o_cvx.x;
   u   = o_cvx.u;
   eta = o_cvx.eta;
   
   x(isnan(x))=0;
   u(isnan(u))=0;
   eta(isnan(eta))=0;
    
%    if o_cvx.status ~= 'Solved'
%        break;
%    end
   
   if o_cvx.delta_norm < 1e-2 && o_cvx.sigma_norm < 1e-2 && o_cvx.nu_norm < 1e-7
       converged = true;
   end
   
   weights.w_dxu = weights.w_dxu*1.5;
   
end

%% plots
close all;

dt = eta/(K-1);
time = (0:K-1) * dt;

figure;
plot(time, x(1,:).*lander_nd.UM)
title('mass')

figure;
plot(time, x(2,:).*lander_nd.UR); hold on;
plot(time, x(3,:).*lander_nd.UR); hold on;
plot(time, x(4,:).*lander_nd.UR); hold on;
title('distance')
legend('x','y','z')

figure;
plot(time, x(5,:).*lander_nd.UV); hold on;
plot(time, x(6,:).*lander_nd.UV); hold on;
plot(time, x(7,:).*lander_nd.UV); hold on;
title('velocity')
legend('x','y','z')

figure;
plot(time, x(11,:)/lander_nd.UT); hold on;
plot(time, x(12,:)/lander_nd.UT); hold on;
plot(time, x(13,:)/lander_nd.UT); hold on;
title('omega rate')
legend('w1','w2','w3')

figure;
plot(time, u(1,:).*lander_nd.UF); hold on;
plot(time, u(2,:).*lander_nd.UF); hold on;
plot(time, u(3,:).*lander_nd.UF); hold on;
plot(time, vecnorm(u).*lander_nd.UF, '--'); hold on;
title('inputs')
legend('x','y','z')

figure;
buffer = 0.1 * max(x(2:4,1));
plot3(x(2,:), x(3,:), x(4,:)); hold on;
plot3(x(2,1), x(3,1), x(4,1),'go'); hold on;
plot3(x(2,end), x(3,end), x(4,end),'ro'); hold on;
grid on;
title(['Trajectory (m) over time']);
xlabel('X distance (m)')
ylabel('Y distance (m)')
zlabel('Z distance (m)')
axis equal;
thelims_max = max([max(x(2,:)) max(x(3,:)) max(x(4,:))]);
thelims_min = min([min(x(2,:)) min(x(3,:)) min(x(4,:))]);
xlim([min(x(2,:))-buffer max(x(2,:))+buffer])
ylim([min(x(3,:))-buffer max(x(3,:))+buffer])
zlim([0 max(x(4,:))+buffer])

