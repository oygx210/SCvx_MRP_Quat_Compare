function [A_bar, B_bar, C_bar, S_bar, z_bar] = discretized_dynamics(vehicle, x, u, sigma)
%discretized_dynamics Generates discretized linear dynamical system
%matrices
%   Does some cool shit, pay attention fools



 
    for k = 0:K-1
        
    end


end

% 	def discretized_dynamics(self, x, u, sigma):
% 		for k in range(self.K -1):
% 			# TODO
% 			self.lin_dyn_0[self.x_ind] = x[:,k];
% 			lin_dyn = np.array(odeint(self.RHS, self.lin_dyn_0, \
% 			(0, self.dt), args=(u[:, k], u[:, k+1], sigma))[1, :]);
% 
% 			B_tau = lin_dyn[self.B_bar_ind].reshape((self.x_n, self.u_n));
% 			C_tau = lin_dyn[self.C_bar_ind].reshape((self.x_n, self.u_n));
% 			S_tau = lin_dyn[self.S_bar_ind];
% 			z_tau = lin_dyn[self.z_bar_ind];
% 
% 
% 			STM = lin_dyn[self.A_bar_ind].reshape([self.x_n, self.x_n]);
% 			self.A_bar[:,k] = STM.flatten(order = 'F');
% 			self.B_bar[:,k] = np.matmul(STM, B_tau).flatten(order = 'F');
% 			self.C_bar[:,k] = np.matmul(STM, C_tau).flatten(order = 'F');
% 			self.S_bar[:,k] = np.matmul(STM, S_tau);
% 			self.z_bar[:,k] = np.matmul(STM, z_tau);
% 		return self.A_bar,self.B_bar,self.C_bar,self.S_bar,self.z_bar;