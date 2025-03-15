function [RK,Quadrature,Model] = getInput()
%% GIVE INPUT PARAMETER
% Sample Input File for Beam Problem 

%% (1) Material
% Linear Elasticity
% Lame Parameters for Young's modulus and poisson ratio
Model.E = 6.1E2; 
Model.nu = 0.49;  % user input
Model.Condition = 'PlaneStress';  % user input: 'PlaneStress' or 'PlaneStrain'
Model.ElasticTensor = getElasticTensor(Model.E,Model.nu,Model.Condition);
Model.DOFu = 2;                  % two dimensional problem

%% (2) Geometry
% rectangular polygon
    
% x1_vertices = [0 48 48 0]';  % user input
% x2_vertices = [-6 -6 6 6]';  % user input
% [x1_vertices,x2_vertices] = interpm([x1_vertices; x1_vertices(1)],...
%                                     [x2_vertices; x2_vertices(1)],2);
% x1_vertices(end) = []; x2_vertices(end) = [];

ym = 2.5;
ymm = 5;
xm = 30;
xmd = 30;
% coff1 = polyfit([0,xmd,xm],[ ym, ymm, ym],2);
% coff2 = polyfit([0,xmd,xm],[-ym,-ymm,-ym],2);
% x1 = linspace(0,xm,8)';
% y1 = coff1(1).*x1.^2 + coff1(2).*x1 + coff1(3);
% x2=flip(x1);
% y2 = coff2(1).*x2.^2 + coff2(2).*x2 + coff2(3);
% 
% % ensure the boundary segments to be counter clockwise
% % [x1_vertices, x2_vertices] = poly2ccw([x1_vertices;x;x], [x2_vertices;y1;y2]);
% [x1_vertices, x2_vertices] = poly2ccw([x1;x2], [y1;y2]);
x1_vertices = [0 xm xm 0]';  % user input
x2_vertices = [-ymm -ymm ymm ymm]';  % user input
[x1_vertices, x2_vertices] = poly2ccw(x1_vertices, x2_vertices);
Model.xVertices = [x1_vertices, x2_vertices];
Model.DomainArea = polyarea(x1_vertices,x2_vertices);

%% (3) Boundary condition
% If an edge is not specified, natural BC with zero traction is imposed.
Model.CriteriaEBC = @(x1,x2) find(x1<=0+1E-7);    % user input
Model.CriteriaNBC = @(x1,x2) find(x1>=xm-(1E-7));   % user input

% beta parameter for Nitches Method
Model.Beta_Nor = 1E2;

% % For verification purpose, provide the exact displacement solution
% syms x1 x2        % use x1 and x2 as x- & y- coordinates
% H = 12; L = 48; D = H; trac = 0.1;   % user input
% I_inertia = (H^3)/12;    % user input
% E = Model.E; nu = Model.nu;

% exact solution
% u1 = (trac*x2/(6*E*I_inertia)).*((6*L-3*x1).*x1+(2+nu)*(x2.^2-(D^2)/4));
% u2 = -(trac/(6*E*I_inertia)).*(3*nu*x2.^2.*(L-x1)+...
%                                (4+5*nu).*((D^2*x1)/4)+(3*L-x1).*x1.^2);
% 
% Model.ExactSolution.u_exact =[u1;u2;];
     
if isfield(Model,'ExactSolution')    %  True if field is in structure array.
    [Model.ExactSolution.S,...
     Model.ExactSolution.g,...
     Model.ExactSolution.t,...
     Model.ExactSolution.b] = getBoundaryConditions(Model);
     Model.ExactSolution.Exist = 1;
else
     Model.ExactSolution.S = @getSebc; 
     Model.ExactSolution.g = @getGebc; 
     Model.ExactSolution.t = @getTraction; 
     Model.ExactSolution.b = @getBodyForce; 
     Model.ExactSolution.Exist = 0;
end

%% (4) Discretization Method
% For general purpose, one can always use A
Model.Discretization.Method = 'A'; 

% (...A) MATLAB built-in FE mesh generator: Default
% if A is chosen, 
% define Hmax: max nodal distance for MATLAB built-in mesh generator
% if Hmax is <=0, then mesh is generated automatically 
% by MATLAB built-in routine
% https://www.mathworks.com/help/pde/ug/pde.pdemodel.generatemesh.html
Model.Discretization.Hmax = 2.4; 


% (...B) Uniform/Non-uniform discretization for rectangular domain:
% if A is chosen, 
% define nx and ny: nx*ny is total number of nodes
% Randomness can be introduced to nodal distribution. 
Model.Discretization.nx1 = 32; 
Model.Discretization.nx2 = 8;
Model.Discretization.Randomness = 0.5;  % 0~1


% (...C) Shestakov distorted discretization for the plate problem:
% if C is chosen, 
% define nc and randomness: 
% nc>1 denotes a refinement parameter
% 0<Distortion<=0.5 is the distortion level parameter.
% This meshing option is specialized for "plate with a hole" geometry
% modifications are needed for other geometries
% as explained in the reference paper 
Model.Discretization.nc = 3;             % >=1
Model.Discretization.Distortion = 0.1;   % 0~0.5


%% (5) RK shape function parameter
RK.KernelFunction = 'SPLIN3';       % SPLIN3
RK.KernelGeometry = 'CIR';          % CIR, REC
RK.NormalizedSupportSize = 2.01;    % suggested order n + 1;
RK.Order = 'Linear';                % Constant, Linear, Quadratic

%% (6) Quadrature rule
Quadrature.Integration = 'SCNI';       % GAUSS, SCNI, DNI
Quadrature.Stabilization = 'N';        % M, N, WO
Quadrature.Option_BCintegration = 'NODAL'; % NODAL OR GAUSS
Quadrature.nGaussPoints = 6; % nGaussPoints per cell
Quadrature.nGaussCells = 5; % nGaussCells on the short side of the domain

%% (7) Plot the variables one prefer
Model.Plot.Discretization = 1;  % plot discretization, nodal representative domain, or Gauss cell
Model.Plot.Displacement = 1;    % plot displacement
Model.Plot.Strain = 1;        % plot strain
Model.Plot.Stress = 1; % plot stress
Model.Plot.DeformedConfiguration = 1;  % plot deformed configuration
Model.Plot.Error = 0; % plot absolute error

end


%% get Elastic Tensor
function [C] = getElasticTensor(E,nu,Condition)
%GETENU 
% Input: Depending on the material type, one may have different input
% Output: Stiffnesss Matrix

switch Condition
    case 'PlaneStrain'
        % Plain Strain
        C =(E/((1+nu)*(1-2*nu)))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, 1-2*nu;];
    case 'PlaneStress'
        % Plain Stress
        C = E / (1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    otherwise
        % Plain Strain
        disp('Default Stiffness Condition: Plane Strain')
        C =(E/((1+nu)*(1-2*nu)))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, 1-2*nu;];
end

end

%% Create Boundary Conditions and body forces
function [ t ] = getTraction(x1,x2,n1,n2)
% Input: 
%   x1 x2: Cartesian coordinate
%   n1,n2: Normal vector at X on boundary
% Output: 
%       t: a 2 by 1 vector for the traction
% ym = 2.5;
% ymm = 5;
% xm = 60;
% xmd = 30;
% coff1 = polyfit([0,xmd,xm],[ ym, ymm, ym],2);
% r = coff1(1).*x1.^2 + coff1(2).*x1 + coff1(3);
% D=r;
H = 15; L = 30; D = H;
I_inertia = H^3/12;
% I_inertia = pi*(2.*r).^4/64;


trac = abs(100);
E = 6.1E1; nu = 0.3;
Condition = 'PlaneStress'; % PlaneStress, or PlaneStrain
C = getElasticTensor(E,nu,Condition);
% epsilon_11 = -(trac*x2.*(6*xm - 6*x1))/(6*E*I_inertia);
% epsilon_22 = (nu*trac*x2.*(xm - x1))/(E*I_inertia);
% epsilon_12 = (trac*x2.^2*(nu + 2))/(3*E*I_inertia);

epsilon_11 = 30*(trac * x1) /(6*E*I_inertia); % 拉伸應變
epsilon_22 = -nu * epsilon_11;     % 泊松效應
epsilon_12 = 0;                     % 沒有剪切應變
% disp(epsilon_11)
stress = C*[epsilon_11;epsilon_22;epsilon_12;];
eta = [n1 0 n2; 0 n2 n1;];
t = eta*stress;
end

function [ SEBC ] = getSebc(x1,x2)
%  Input: 
%  x1 x2: Cartesian coordinate
% Output: 
%   SEBC: a 2 by 2 matrix for switch on EBC
SEBC = diag([1 1]); 
end

function [ gEBC ] = getGebc(x1,x2)
%  Input: 
%  x1 x2: Cartesian coordinate
% Output: 
%   gEBC: a 2 by 1 vector of prescribed displacement on EBC
ym = 2.5;
ymm = 15;
xm = 30;
xmd = 30;
coff1 = polyfit([0,xmd,xm],[ ym, ymm, ym],2);
r = coff1(1).*x1.^2 + coff1(2).*x1 + coff1(3);
I = pi*(2.*r).^4/64;

trac = abs(1000);
E = 6.1E2; nu = 0.49;

% u1_exact = (trac*x2/(6*E*I)).*((6*L-3*x1).*x1+(2+nu)*(x2.^2-(H^2)/4));
% u2_exact = -(trac/(6*E*I)).*(3*nu*x2.^2.*(L-x1)+(4+5*nu).*((H^2*x1)/4)+(3*L-x1).*x1.^2);
% gEBC = [u1_exact;
%          u2_exact];
Condition = 'PlaneStress'; % PlaneStress, or PlaneStrain
C = getElasticTensor(E,nu,Condition);
epsilon_11 = (trac * x1) /(6*E*I); % 拉伸應變
epsilon_22 = -nu * epsilon_11;     % 泊松效應
epsilon_12 = 0;                     % 沒有剪切應變
% disp(epsilon_11)
stress = C*[epsilon_11;epsilon_22;epsilon_12;];
gEBC = [stress(1);stress(2)];
end

function [ b ] = getBodyForce(x1,x2)
% Input: 
%  x1 x2: Cartesian coordinate
% Output: 
%      b: a 2 by 1 vector for the body force
b = [0; -0.08;];
end