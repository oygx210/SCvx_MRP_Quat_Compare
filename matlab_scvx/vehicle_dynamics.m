classdef vehicle_dynamics
    %vehicle_dynamics General class for vehicle dynamics and discretization
    %   This class will allow you to symbollically create input output
    %   functions for jacobians, and then create their instances.
    
    properties
        vehicle
        f_sym, f_handle
        A_sym, A_handle
        B_sym, B_handle
        x_sym
        u_sym
        K, dt, m, n
        x_ind
        A_bar_ind, A_bar
        B_bar_ind, B_bar
        C_bar_ind, C_bar
        S_bar_ind, S_bar
        z_bar_ind, z_bar
        lds_0
    end
    
    methods
        function self = vehicle_dynamics(lander)
            %vehicle_dynamics Construct an instance of the vehicle dynamics
            %class
            %   Detailed explanation goes here
            self.vehicle = lander;
            [f, A, B, x, u] = generate_vehicle_jacobians(lander);
            self.f_sym = f;
            self.f_handle = matlabFunction(f,'Vars', {[x.' u.'].'});
            self.A_sym = A;
            self.A_handle = matlabFunction(A,'Vars', {[x.' u.'].'});
            self.B_sym = B;
            self.B_handle = matlabFunction(B,'Vars', {[x.' u.'].'});
            self.x_sym = x;
            self.u_sym = u;
            self.K = lander.K;
            self.dt = 1/(lander.K-1);
            self.m = lander.m;
            self.n = lander.n;
            
            self.x_ind = 1:self.m;
            self.A_bar_ind = self.m + (1:(self.m^2));
            self.B_bar_ind = self.A_bar_ind(end) + (1 : (self.m*self.n));
            self.C_bar_ind = self.B_bar_ind(end) + (1 : (self.m*self.n));
            self.S_bar_ind = self.C_bar_ind(end) + (1 : self.m);
            self.z_bar_ind = self.S_bar_ind(end) + (1 : self.m);
            
            self.A_bar = zeros([self.m self.m self.K-1]);
            self.B_bar = zeros([self.m self.n self.K-1]);
            self.C_bar = zeros([self.m self.n self.K-1]);
            self.S_bar = zeros([self.m self.K-1]);
            self.z_bar = zeros([self.m self.K-1]);
            
            self.lds_0 = zeros(self.z_bar_ind(end), 1);
            self.lds_0(self.A_bar_ind) = reshape(eye(self.m),[],1);
        end
        
        % Compute an instance of A, B, and the nonlinear dynamics
        function [f,A,B] = compute_matrices(self, x, u)
            f = self.f_handle([x; u]).';
            A = self.A_handle([x; u]);
            B = self.B_handle([x; u]);
        end
        
        % Generate the RHS of the system matrices that we'll need to
        % integrate later. We should only integrate from 0 to dtau
        % Dynamics are time invariant.
        function dLdt = RHS(self, t, lds, u_k, u_kp, sigma)
            dLdt = zeros(size(lds));
            ak = (self.dt - t)/ self.dt;
            bk = t/self.dt;
            u = (ak*u_k) + (bk*u_kp);
            
            % Extract the state information
            x = lds(self.x_ind);
            % Pull out the previous intermediate state transition matrix
            Phi_xi_tk = reshape(lds(self.A_bar_ind), [self.m, self.m]);

            [f_xi, A, B] = compute_matrices(self, x, u);
            A_xi = sigma*A;
            B_xi = sigma*B;
            f_xi(isnan(f_xi))=0;
            A_xi(isnan(A_xi))=0;
            B_xi(isnan(B_xi))=0;
            
            z_xi = -A_xi*x - B_xi*u;
            B_mult = reshape(Phi_xi_tk\B_xi, [], 1);
            
            dLdt(self.x_ind)    = (sigma * f_xi).';
            dLdt(self.A_bar_ind)= reshape(A_xi*Phi_xi_tk, [], 1);
            dLdt(self.B_bar_ind)= reshape(ak*B_mult, [], 1);
            dLdt(self.C_bar_ind)= reshape(bk*B_mult, [], 1);
            dLdt(self.S_bar_ind)= (Phi_xi_tk\f_xi).';
            dLdt(self.z_bar_ind)= (Phi_xi_tk\z_xi).';
        end
        
        function output_struct = discretized_dynamics(self, x, u, sigma)
            odesettings = odeset('AbsTol', 1E-12, 'RelTol', 1E-12);
            for k = 1:self.K-1
                % State history over all time (mx1 each step)
                self.lds_0(self.x_ind) = x(:,k);
                nldynamics = @(t, lds) self.RHS(t, lds, u(:,k), u(:,k+1), sigma);
%               nldynamics = @(t, lds) self.RHS(t, self.lds_0, u(:,k), u(:,k+1), sigma);
                lds_out = ode45(nldynamics, [0 self.dt], self.lds_0,odesettings);
                lds_out = lds_out.y(:, end);
                
                % Pull ou the values and reshape them
                A_tau = reshape(lds_out(self.A_bar_ind), [self.m self.m]);
                B_tau = reshape(lds_out(self.B_bar_ind), [self.m self.n]);
                C_tau = reshape(lds_out(self.C_bar_ind), [self.m self.n]);
                S_tau = lds_out(self.S_bar_ind);
                z_tau = lds_out(self.z_bar_ind);
                
                A_tau(isnan(A_tau))=0;
                B_tau(isnan(B_tau))=0;
                C_tau(isnan(C_tau))=0;
                S_tau(isnan(S_tau))=0;
                z_tau(isnan(z_tau))=0;
                
                % Fill in the next STM
                self.A_bar(:,:,k) = A_tau;
                self.B_bar(:,:,k) = A_tau*B_tau;
                self.C_bar(:,:,k) = A_tau*C_tau;
                self.S_bar(:,k)   = A_tau*S_tau;
                self.z_bar(:,k)   = A_tau*z_tau;
            end
            output_struct.A_bar = self.A_bar;
            output_struct.B_bar = self.B_bar;
            output_struct.C_bar = self.C_bar;
            output_struct.S_bar = self.S_bar;
            output_struct.z_bar = self.z_bar;
        end
    end
end





















