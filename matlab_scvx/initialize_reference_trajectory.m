function [x_0, u_0] = initialize_reference_trajectory(lander, K)
%initialize_reference_trajectory Generates initial solve trajectory
%   Outputs initial solution for trajectory and control inputs over K

    X0 = lander.X0;
    U0 = lander.U0;
    XT = lander.XT;
    
    x_0 = zeros(length(X0), K);
    u_0 = zeros(length(U0), K);
    
    % all of this shit is naive frankly.
    Km = K - 1;
    for k = 0:Km
        a1 = (Km  - k)/Km;
        a2 = k/Km;
        m_k = X0(1);
        r_k = a1*X0(2:4) + a2*XT(2:4);
        v_k = a1*X0(5:7) + a2*XT(5:7);
        s_k = X0(8:10);
        w_k = X0(11:13);
        x_0(:,k+1) = [m_k; r_k; v_k; s_k; w_k];
        % hover the beast
        u_k = -lander.g .*m_k;
        u_0(:,k+1) = [0 0 min(max(u_k(3), lander.Fth_min), lander.Fth_max)].';
    end
end

