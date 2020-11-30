function output_cvx = scvx_subproblem(vehicle, lds, x_last, u_last, eta_last, weights)
%scvx_subproblem Solves one iteration of a SCVx problem
%   This requires knowledge of the dynamics and previousl control and input
%   histories. It returns a newly optimized trajectory and input control
%   history with free final time solution and other weights, timing, and
%   error statistics of that particular solution.
%   Note that Sigma was replaced with Eta due to symbolic conflict in cvx.

% Break things out to make it cleaner.
K = vehicle.K;
m = vehicle.m;
n = vehicle.n;
A_bar = lds.A_bar;
B_bar = lds.B_bar;
C_bar = lds.C_bar;
S_bar = lds.S_bar;
z_bar = lds.z_bar;
w_nu  = weights.w_nu;
w_dxu = weights.w_dxu;
w_ds  = weights.w_ds;
w_s   = weights.w_s;

%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN THE BEGUINE
cvx_solver SDPT3
cvx_begin quiet
    % Define our variables
    variable x(m,K)
    variable u(n,K)
    variable eta nonnegative
    variable nu(m,K-1)
    variable delta_norm nonnegative
    variable eta_norm nonnegative

    % Objective Function
    minimize( w_s*eta + w_nu*norm(nu, Inf) + w_dxu*delta_norm + w_ds*eta_norm )
%     minimize( -x(1,K) + w_nu*norm(nu, Inf) + w_dxu*delta_norm + w_ds*eta_norm )

    % Constraints and Dynamics
    subject to
        % Constrain the initial and terminal conditions. Final Mass Free.
        x(:,1)      == vehicle.X0;
        x(2:end, K) == vehicle.XT(2:end);
        x(1,:) >= vehicle.m_dry;

        % Dynamics
        for  k = 1:K-1
            x(:,k+1) == A_bar(:,:,k)*x(:,k) + ...
                        B_bar(:,:,k)*u(:,k) + ...
                        C_bar(:,:,k)*u(:,k+1) + ...
                        S_bar(:,k).*eta + ...
                        z_bar(:,k) + ...
                        nu(:,k);
            % Mass should Never Increase
            x(1,k+1) <= x(1,k);
        end
        
        % Positive time
        eta >= 0.0001;
        
        % Trust Region Constraints
        delta_x = sum(x - x_last)*sum(x - x_last).';
        delta_u = sum(u - u_last)*sum(u - u_last).';
        delta_x + delta_u    <= delta_norm;

        delta_eta = (eta - eta_last);
        norm(delta_eta, Inf) <= eta_norm;
        
        % Original Way of doing these constraints
        % delta_x = sum_square(x - x_last);
        % delta_u = sum_square(u - u_last);
        % norm(delta_x + delta_u, 1) <= delta_norm;
        
        
        % Glideslope, Max Angle, Max Omega, Gimbal Angle, Upper Thrust,
        % Linearized Lower Thrust.
        %    1  2  3  4  5  6  7  8  9 10 11 12 13
        %   [m r1 r2 r3 v1 v2 v3 s1 s2 s3 w1 w2 w3]
        for  k = 1:K
            % Check the order of NED here
            norm(x(2:3, k))   <= x(4,k)/vehicle.glideslope;
%             norm(x(11:13,k))  <= vehicle.omega_max;
            norm(x(8:10, k))  <= vehicle.tan_theta_max;
            cos(vehicle.gimbal_max)*norm(u(:,k))  <= u(3, k);
            u_min_lin = (u_last(:, k)./norm(u_last(:, k))).' * u(:, k);
            u_min_lin         >= vehicle.Fth_min;
			norm(u(:,k))      <= vehicle.Fth_max;
        end

cvx_end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_cvx.x = x;
output_cvx.u = u;
output_cvx.eta = eta;
output_cvx.delta_norm = delta_norm;
output_cvx.sigma_norm = eta_norm;
output_cvx.nu_norm = norm(nu, Inf);
output_cvx.status = cvx_status;

end

