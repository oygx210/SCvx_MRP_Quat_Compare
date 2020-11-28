function [f,A,B,x,u] = generate_vehicle_jacobians(vehicle)
%generate_vehicle_jacobians Summary of this function goes here
%   Detailed explanation goes here

    
    syms m r1 r2 r3 v1 v2 v3 s1 s2 s3 w1 w2 w3 u1 u2 u3 real    
    f = sym('f',[1 13]);
    %    1  2  3  4  5  6  7  8  9 10 11 12 13
    x = [m r1 r2 r3 v1 v2 v3 s1 s2 s3 w1 w2 w3].';
    u = [u1 u2 u3].';

    sigma = x(8:10);
    I_c   = vehicle.I_c;
    r_com = vehicle.r_com;
    C_N2B = MRP2C(sigma).';
    g_planet_N = vehicle.g;
    
    % Nonlinear dynamics (constant gravity, not Newtonian)
    % mdot v1 v2 v3 a1 a2 a3 sd1 sd2 sd3 wd1 wd2 wd3
    f(1) = -vehicle.alpha_mdot * norm(u);
    f(2:4) = x(5:7);
    f(5:7) = ((C_N2B * u)/x(1)) + g_planet_N;
    mrp_B = ((1- sigma.'*sigma)*eye(3) + 2*skew(sigma) + 2*(sigma*sigma.'));
    f(8:10) = (1/4).*mrp_B*x(11:13);
    f(11:13)= I_c \ (skew(r_com)*u - skew(x(11:13))*I_c*x(11:13));
    
    f = simplify(f);
    A = simplify(jacobian(f,x));
    B = simplify(jacobian(f,u));
end

