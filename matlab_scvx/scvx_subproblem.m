function [x_cvx, u_cvx, sigma_cvx, delta_norm_cvx, ...
        sigma_norm_cvx, nu_cvx, error, stats] = scvx_subproblem(lander, A_bar, ...
                                        B_bar, C_bar, Sigma_bar, z_bar, ...
                                        x, u, sigma, weights, K);
%scvx_subproblem Solves one iteration of a SCVx problem
%   This requires knowledge of the dynamics and previousl control and input
%   histories. It returns a newly optimized trajectory and input control
%   history with free final time solution and other weights, timing, and
%   error statistics of that particular solution.

    cvx_solver ECOS
    cvx_begin
        
        % Define our variables
        variables u(3,N) z(1,N) s(1,N) r(3,N) v(3,N)
        
        % Objective Function
        minimize(-z(N))			% objective function
        
        % Constraints and Dynamics
        subject to
            r(:,1) == r_0;
            v(:,1) == v_0;
            r(:,N) == r_N;
            v(:,N) == v_N;
            z(1) == log(m_t);

            for  k = 1:N-1
                r(:,k+1) == r(:,k) + ((dt/2)*(v(:,k) + v(:,k+1))) +(((dt^2)/12)*(u(:,k+1) - u(:,k)));
            end

            for k=1:N
                norm(u(:,k)) <= s(1,k);
            end

    cvx_end

end

