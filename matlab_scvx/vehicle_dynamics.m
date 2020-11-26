classdef vehicle_dynamics
    %vehicle_dynamics General class for vehicle dynamics and discretization
    %   This class will allow you to symbollically create input output
    %   functions for jacobians, and then create their instances.
    
    properties
        vehicle
        f_sym
        A_sym
        B_sym
        x_sym
        u_sym
        K
        dt
        m
        n
        x_ind
        A_bar_ind
        B_bar_ind
        C_bar_ind
        S_bar_ind
        z_bar_ind
    end
    
    methods
        function self = vehicle_dynamics(lander)
            %vehicle_dynamics Construct an instance of the vehicle dynamics
            %class
            %   Detailed explanation goes here
            self.vehicle = lander;
            [f, A, B, x, u] = generate_vehicle_jacobians(lander);
            self.f_sym = f;
            self.A_sym = A;
            self.B_sym = B;
            self.x_sym = x;
            self.u_sym = u;
            self.K = lander.K;
            self.dt = 1/(lander.K-1);
            self.m = lander.m;
            self.n = lander.m;
            
            self.x_ind = 1:self.m;
            self.A_bar_ind = self.m + (1:(self.m^2));
            self.B_bar_ind = self.A_bar_ind(end) + (1 : (self.m*self.n));
            self.C_bar_ind = self.B_bar_ind(end) + (1 : (self.m*self.n));
            self.S_bar_ind = self.C_bar_ind(end) + (1 : self.m);
            self.z_bar_ind = self.S_bar_ind(end) + (1 : self.m);
        end
        
        % Compute an instance of A, B, and the nonlinear dynamics
        function [f,A,B] = compute_matrices(self, x, u)
            sym_vars = [self.x_sym',self.u_sym']';
            numerical_values = [x',u']';
            f = double(subs(self.f_sym, sym_vars, numerical_values));
            A = double(subs(self.A_sym, sym_vars, numerical_values));
            B = double(subs(self.B_sym, sym_vars, numerical_values));
        end
        
        % Generate the RHS of the system matrices that we'll need to
        % integrate later. We should only integrate from 0 to dtau
        % Dynamics are time invariant.
        function dLdt = RHS(self, lds, t, u_k, u_kp, sigma)
            ak = (self.dt - t)/ self.dt;
            bk = t/self.dt;
            u = (ak*u_k) + (bk*u_kp);
            
            x = lds(self.x_ind);
            Phi_k_km = reshape(lds(self.A_bar_ind), [self.m, self.m]);


            [f_xi,A,B] = compute_matrices(self, x, u);
            
            A_xi = sigma*A;
            B_xi = sigma*B;
            
            f_xi(isnan(f_xi))=0;
            A_xi(isnan(A_xi))=0;
            B_xi(isnan(B_xi))=0;
            
            
            
            dLdt(self.x_ind) = (sigma * f_xi).';
            dLdt(self.A_bar_ind)= reshape(A_xi*last_Phi     , [], 1);
            dLdt(self.B_bar_ind)= reshape(ak.*Phi_tk_xi*B_xi, [], 1);
            dLdt(self.C_bar_ind)= reshape(bk.*Phi_tk_xi*B_xi, [], 1);
            dLdt(self.S_bar_ind)= (Phi_tk_xi*f_xi).';
            dLdt(self.z_bar_ind)= reshape( Phi_tk_xi*z_xi , [], 1);
        end
        
        function [A_bar,B_bar,C_bar,S_bar,z_bar] = discretized_dynamics(self, x, u, sigma)
            
            for k = 1:self.K
               A_bar = zeros(self.m, self.m);
               B_bar = zeros(self.m, self.n);
               C_bar = zeros(self.m, self.n);
               S_bar = zeros(self.m, 1);
               z_bar = zeros(self.m, 1);
            end
        end

        
        
    end
end





















