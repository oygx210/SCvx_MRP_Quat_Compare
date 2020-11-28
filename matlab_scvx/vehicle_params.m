% Vehicle configuration

lander.m = 300;
lander.g = [0 0 -9.81].';            
lander.Fth_max = 3113.76;
lander.Fth_min = 0.15 * lander.Fth_max;
lander.Isp = 180; % purposely unrealistic
lander.alpha_mdot = 1/(norm(lander.g)*lander.Isp);
lander.m_wet = lander.m;
lander.m_dry = 50;

h = 2.5;
r = 0.6;
% xyz with z going through the boresight vector
lander.I_c = diag([(1/12)*lander.m*(3*r*r + h*h), ...
			(1/12)*lander.m*(3*r*r + h*h) (1/2)*lander.m*r*r]);

lander.r_com = [0 0 -1.2191].';

% constraints
lander.glideslope = tan(deg2rad(20));
lander.gimbal_max = deg2rad(15);
lander.cos_theta_max = cos(deg2rad(89));
lander.tan_theta_max = tan(deg2rad(89)/4);
lander.omega_max = deg2rad(10);
