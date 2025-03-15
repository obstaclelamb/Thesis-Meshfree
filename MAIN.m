
clear all;close all;close all hidden;close all force; format short e;clc;

% add path to subroutine folder 
path_MAINCODE = 'MAIN_PROGRAM'; 
rmpath(path_MAINCODE); addpath(path_MAINCODE);

% add path to input file folder
% path_problem = '01_LinearPatchTest'; 
% path_problem = '02_BeamProblem'; 
% path_problem = '03_PlateProblem'; 
% % path_problem = '04_PoissonProblem_src'; % is the diffusion problem in paper; 
% path_problem = '05_FEA_Model'; 
% path_problem = '06_TensileTest'; 
path_problem = 'X_vary_BeamProble'; % test for muscle
rmpath(path_problem); addpath(path_problem)


%% (1) Pre-Process
% (1.1) Generate input files (RK parameter,Quadrature, Material, BC)
[RK,Quadrature,Model] = getInput();

% (1.2) Generate nodal discretization (three options available)
[Discretization,Model] = Pre_GenerateDiscretization(Model);   

% (1.3) Evaluate support size & output model set-up parameters
[RK,Quadrature] = Pre_Initialization(RK,Quadrature,Discretization);

% (1.4) Generate quadrature weight & evaluatioin point
[Quadrature] = Pre_GenerateQuadraturePoints(Quadrature,Discretization);

% (1.5) Evaluate RK shape functions
[Quadrature] = Pre_GenerateShapeFunction(RK,Quadrature);

%% (2) Main Solver
% (2.1) Assemble global stiffness matrix & force vector
[K,F] = MatrixAssembly(Quadrature,Model);

% (2.2) Solve linear algebraic equation
RK.dI = mldivide(K,F);  % MATLAB built-in matrix solver i.e., dI = K\F

%% (3) Post-Process
% (3.1) visualization
[uhI,Strain,Stress,eru,erdu] = PostProcess(RK,Quadrature,Model,Discretization);

% clean path
rmpath(path_MAINCODE)
rmpath(path_problem)