function [RK,Quadrature] = Pre_Initialization(RK,Quadrature,Discretization)
% This subroutine initializes the matrices and vectors for the simulation. 
% It also displays the input information in the MATAB command window

%% Read coordinates from discretization and save it to RK structure
xI = Discretization.xI; nP = Discretization.nP;
RK.xI = Discretization.xI; RK.nP = Discretization.nP;

%% Display the parameter choose
disp(['Start Simulation: the domain is partition to: ',num2str(Discretization.nP),' RK nodes'])

%% Display the parameter choose for order
% if ~exist('RK.Order');  RK.Order = 'Linear'; end;
switch RK.Order
    case {'Constant'}
        disp(['Order of Basis: ',RK.Order])
    case {'Linear'}
        disp(['Order of Basis: ',RK.Order])
    case {'Quadratic'}    
        disp(['Order of Basis: ',RK.Order])
        disp('Note: SCNI cannot pass quadratic patch test.')
    otherwise
        RK.Order = 'Linear';
        disp('Your choice of order is inconsistent. Swtich to defualt basis.')
        disp(['Defualt order of basis: ',RK.Order])
end

%% display kernel
% if ~exist('RK.KernelFunction');  RK.KernelFunction = 'CBSPLIN'; end;
% if ~exist('RK.NormalizedSupportSize');  RK.NormalizedSupportSize = 2.001; end;

disp(['Kernel Function Type: ',RK.KernelFunction])
disp(['Normalized Support Size: ',num2str(RK.NormalizedSupportSize)])

%% Base on the kernel geometry (rectangular or circular )
% if ~exist('RK.KernelGeometry');  RK.KernelGeometry = 'CIR'; end;
switch RK.KernelGeometry
    case 'CIR'
        SupportDistanceType = 'euclidean';
        disp('Support Type: Circular')
    case 'REC'
        disp('Support Type: Rectangular')
        SupportDistanceType = 'chebychev';    
    otherwise
        RK.KernelGeometry = 'CIR'; 
        disp('Default Support Type: Circular')
        SupportDistanceType = 'euclidean';
end
RK.SupportDistanceType = SupportDistanceType;

%% Display the parameter chosen for stabilization technique
% if ~exist('Quadrature.Integration');  Quadrature.Integration = 'SCNI'; end;
% if ~exist('Quadrature.Option_BCintegration'); Quadrature.Option_BCintegration = 'NODAL'; end;
switch Quadrature.Integration
    case {'SCNI','DNI'}
        disp(['Integration Method: ',Quadrature.Integration])
        disp(['Boundary Integration Method: ',Quadrature.Option_BCintegration])
    case 'GAUSS'
        Quadrature.Option_BCintegration = 'GAUSS';
        disp(['Integration Method: ',Quadrature.Integration])  
        disp(['Background Integration Cell: ',num2str(Quadrature.nGaussPoints),' by ',num2str(Quadrature.nGaussPoints)])  
    otherwise
        Quadrature.Integration = 'SCNI';
        Quadrature.Option_BCintegration = 'NODAL';
        disp(['Defualt Integration Method: ',Quadrature.Integration])
        disp(['Default Boundary Integration Method: ',Quadrature.Option_BCintegration])
end

% if ~exist('Quadrature.Stabilization','var');  Quadrature.Stabilization = 'N'; end;
switch Quadrature.Integration
    case {'SCNI','DNI'}
        switch Quadrature.Stabilization
            case {'MSCNI','MDNI','M'}
                Quadrature.MtypeStablization.SHPDX1_Is = cell(nP,1);
                Quadrature.MtypeStablization.SHPDX2_Is = cell(nP,1);
                Quadrature.MtypeStablization.AREA_Is = cell(nP,1);                        
                disp(['Stabilization Technique: ',Quadrature.Stabilization])
            case {'NSNI','NSCNI','NDNI','N'}
                Quadrature.NtypeStablization.ML = zeros(nP,2);
                M1L = zeros(nP,1);
                M2L = zeros(nP,1);
                Quadrature.NtypeStablization.ML(:,1) = M1L;
                Quadrature.NtypeStablization.ML(:,2) = M2L; 
                Quadrature.NtypeStablization.SHPDX1X1 = sparse(nP,nP); % Two Deravitive in x1
                Quadrature.NtypeStablization.SHPDX2X2 = sparse(nP,nP); % Two Deravitive in x2
                Quadrature.NtypeStablization.SHPDX1X2 = sparse(nP,nP); % Deravitive in x1 and Deravitive in x2
                Quadrature.NtypeStablization.SHPDX2X1 = sparse(nP,nP); % Deravitive in x2 and Deravitive in x1
                disp(['Stabilization Technique: ',Quadrature.Stabilization])
            otherwise
                disp('No stabilization technique is employed.')
        end
    otherwise
        disp('Gauss integration without stabilization technique.')
end


%% Base on the kernel type, we can decide the the k nearest neighbor 
% such that the nodal distance can be determined
% use 5 to ensure at least 4 nearby neighbors for each node 
[~,h_4point] = knnsearch(xI,xI,'K',5,'Distance',SupportDistanceType);  
NodalSpacing = h_4point(:,5); % obtain the largest nodal spacing size, h

% in Nitches method, one need to evaluate the Nitches parameter by spacing,
Quadrature.NodalSpacing = NodalSpacing;

Diala_Support = RK.NormalizedSupportSize*NodalSpacing; 

RK.SupportSize = Diala_Support; 


end