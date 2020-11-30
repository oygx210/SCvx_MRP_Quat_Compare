function vehicle_nd = compute_nd_factors(vehicle, tf_guess, r0, v0, sigma, omega)
%compute_nd_factors Computes nondimensionalized factors for vehicle object
%   Sigma can be whatever attitude formalism you choose.

    vehicle_nd.UM = vehicle.m;
    vehicle_nd.UR = norm(r0);
    vehicle_nd.UT = tf_guess;
    vehicle_nd.UV = vehicle_nd.UR/vehicle_nd.UT;
    vehicle_nd.UA = vehicle_nd.UV/vehicle_nd.UT;
    vehicle_nd.UF = vehicle_nd.UA*vehicle_nd.UM;
    vehicle_nd.UI = vehicle_nd.UM*(vehicle_nd.UR^2);
    
%     vehicle_nd.UM = 1;
%     vehicle_nd.UR = 1;
%     vehicle_nd.UT = 1;
%     vehicle_nd.UV = 1;
%     vehicle_nd.UA = 1;
%     vehicle_nd.UF = 1;
%     vehicle_nd.UI = 1;
    
    vehicle_nd.m = vehicle.m                        / vehicle_nd.UM;
    vehicle_nd.g = vehicle.g                        / vehicle_nd.UA;
    vehicle_nd.Fth_min = vehicle.Fth_min            / vehicle_nd.UF;
    vehicle_nd.Fth_max = vehicle.Fth_max            / vehicle_nd.UF;
    vehicle_nd.Isp = vehicle.Isp                    / vehicle_nd.UT;
    vehicle_nd.alpha_mdot = vehicle.alpha_mdot      * (vehicle_nd.UA*vehicle_nd.UT);
    vehicle_nd.m_wet = vehicle.m_wet                / vehicle_nd.UM;
    vehicle_nd.m_dry = vehicle.m_dry                / vehicle_nd.UM;
    vehicle_nd.I_c = vehicle.I_c                    / vehicle_nd.UI;
    vehicle_nd.r_com = vehicle.r_com                / vehicle_nd.UR;
    
    vehicle_nd.glideslope    = vehicle.glideslope;
    vehicle_nd.gimbal_max    = vehicle.gimbal_max;
    vehicle_nd.cos_theta_max = vehicle.cos_theta_max;
    vehicle_nd.tan_theta_max = vehicle.tan_theta_max;
    vehicle_nd.omega_max     = vehicle.omega_max    *vehicle_nd.UT;
    
    % initial condition for state and input
    vehicle_nd.X0 = [vehicle_nd.m; r0./vehicle_nd.UR; v0./vehicle_nd.UV; ...
                        sigma; omega*vehicle_nd.UT];
                    
    vehicle_nd.U0 = [0 0 0].';

end

