function [function_S,function_g,function_traction,function_b] = getBoundaryConditions(Model)
% This subroutine outputs the functional of the prescribed boundary conditions 
% t=n^TCepsilon^exact, g=u^exact, S=I_(2×2) and the body force b=sigma_(ij,j)^exact 
% according to given expressions of displacement 
% u_exact in terms of x_1 and x_2 in a symbolic format

syms x1 x2 n1 n2

% function handle for essential boundary condition g
u = Model.ExactSolution.u_exact; 
C = Model.ElasticTensor;
function_g = matlabFunction(u);

% function handle for stress
epsilon_x1 = diff(u(1),x1);
epsilon_x2 = diff(u(2),x2);
epsilon_x12 =  (diff(u(1),x2)+diff(u(2),x1));
stress = C*[epsilon_x1;epsilon_x2;epsilon_x12];
% function_stress = matlabFunction(stress,'Vars',[x1 x2]);
% function handle for traction
eta = [n1 0 n2; 0 n2 n1;];
traction = eta*stress;
function_traction = matlabFunction(traction,'Vars',[x1 x2 n1 n2]);

% function handle for body force b
b = [divergence([stress(1),stress(3)],[x1,x2]);
     divergence([stress(3),stress(2)],[x1,x2]);];
function_b = matlabFunction(b,'Vars',[x1 x2]);

% function handle for switch S
function_S = matlabFunction(sym(diag([1 1])),'Vars',[x1 x2]); 


end

