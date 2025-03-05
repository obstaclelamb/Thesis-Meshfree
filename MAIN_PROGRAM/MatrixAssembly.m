function [K_total,F_total] = MatrixAssembly(Quadrature,Model)
% This subroutine assembles the stiffness matrix and force vector
% number of nodal from Quadrature
nP = Quadrature.nP;
% from QDRTURE to get evaluation points for domain
xQuad = Quadrature.Domain.xQuad; 
nQuad = Quadrature.Domain.nQuad;
% load E, nu, C from Model 
E = Model.E; nu = Model.nu;
C = Model.ElasticTensor;
% exact solution function handle
f_S = Model.ExactSolution.S;
f_g = Model.ExactSolution.g;
f_t = Model.ExactSolution.t;
f_b = Model.ExactSolution.b;
% Degrees of freedom
DOFu = Model.DOFu;
% Penalty Parameters for Nitche's Method
% beta = Model.Beta_Nor*E/sqrt((Model.DomainArea/nP)); % Penalty Number
beta = Model.Beta_Nor*E/max(Quadrature.NodalSpacing); % use the max h
% Matrix/Force for tangent stiffness
KIJ_c = sparse ( nP*DOFu , nP*DOFu ) ;
FI_b = sparse( nP*DOFu,1) ; 
FI_t = sparse( nP*DOFu,1) ; 
% Matrix/Force for boundary condition by Nitches method
KIJ_g = sparse ( nP*DOFu , nP*DOFu ) ; 
KIJ_beta = sparse ( nP*DOFu , nP*DOFu ) ;
FI_g = sparse( nP*DOFu,1) ; 
FI_beta = sparse( nP*DOFu,1) ; 
% B matrix and Psi matrix 
B = zeros(3,nP*DOFu) ; 
PSI = zeros(2,nP*DOFu) ; 

%% Domain Integration
Weight = Quadrature.Domain.Weight;
SHP = Quadrature.SHP;
SHPDX1 = Quadrature.SHPDX1;
SHPDX2 = Quadrature.SHPDX2;

% if stabilization is used, load the corresponding stabilization terms
switch Quadrature.Integration
    case {'SCNI','DNI'}
        switch Quadrature.Stabilization
            case {'M'}
            SHPDX1_Is = Quadrature.MtypeStablization.SHPDX1_Is ;
            SHPDX2_Is = Quadrature.MtypeStablization.SHPDX2_Is;
            Area_Is = Quadrature.MtypeStablization.AREA_Is;
            B_Is = zeros(3,nP*DOFu);
            case {'N'}
            SHPDX1X1 = Quadrature.NtypeStablization.SHPDX1X1 ;
            SHPDX2X2 = Quadrature.NtypeStablization.SHPDX2X2;
            SHPDX1X2 = Quadrature.NtypeStablization.SHPDX1X2;
            SHPDX2X1 = Quadrature.NtypeStablization.SHPDX2X1;
            M1N = Quadrature.NtypeStablization.M(:,1);
            M2N = Quadrature.NtypeStablization.M(:,2);
            B_1I_ig = zeros(3,nP*DOFu);
            B_2I_ig = zeros(3,nP*DOFu);
        end
end
% This is the number% bar to show the progress
sub_textprogressbar('Loop Over Quadrature Points, Progress: ')
% loop over integration points
for idx_nQuad = 1:nQuad
    % Progress Bar
    sub_textprogressbar(100*idx_nQuad/nQuad)
    B(:,:) = 0; PSI(:,:) = 0;
    % evalutate B and Psi matrix by shape function
    B(1,1:2:end) = SHPDX1(idx_nQuad,:);
    B(2,2:2:end) = SHPDX2(idx_nQuad,:);
    B(3,1:2:end) = SHPDX2(idx_nQuad,:);
    B(3,2:2:end) = SHPDX1(idx_nQuad,:);
    PSI(1,1:2:end) = SHP(idx_nQuad,:); PSI(2,2:2:end) = SHP(idx_nQuad,:);
    % sitffness matrix assembly
    KIJ_c = KIJ_c + (B')*(C)*B*(Weight(idx_nQuad,1));
    % stabilization terms
    switch Quadrature.Integration
        case {'SCNI','DNI'}
            switch Quadrature.Stabilization
                case {'M'} 
                [J_MSCNI,~] = size(SHPDX1_Is{idx_nQuad}(:,:));
                for j = 1:J_MSCNI % each edge for each cell
                    B_Is(:) = 0;
                    B_Is(1,1:2:end) = SHPDX1_Is{idx_nQuad}(j,:);
                    B_Is(2,2:2:end) = SHPDX2_Is{idx_nQuad}(j,:);
                    B_Is(3,1:2:end) = SHPDX2_Is{idx_nQuad}(j,:);
                    B_Is(3,2:2:end) = SHPDX1_Is{idx_nQuad}(j,:);

                   KIJ_c = KIJ_c + (B_Is'-B')*(C)*(B_Is-B)*Area_Is{idx_nQuad}(j,1);
                end
                case {'N'}
                    B_1I_ig(:) = 0;
                    B_2I_ig(:) = 0;
                    B_1I_ig(1,1:2:end) = (SHPDX1X1(idx_nQuad,:));
                    B_1I_ig(2,2:2:end) = (SHPDX2X1(idx_nQuad,:));
                    B_1I_ig(3,1:2:end) = (SHPDX2X1(idx_nQuad,:));
                    B_1I_ig(3,2:2:end) = (SHPDX1X1(idx_nQuad,:));
                    B_2I_ig(1,1:2:end) = (SHPDX1X2(idx_nQuad,:));
                    B_2I_ig(2,2:2:end) = (SHPDX2X2(idx_nQuad,:));
                    B_2I_ig(3,1:2:end) = (SHPDX2X2(idx_nQuad,:));
                    B_2I_ig(3,2:2:end) = (SHPDX1X2(idx_nQuad,:));
                    KIJ_c = KIJ_c + ((B_1I_ig)'*(C)*(B_1I_ig)*M1N(idx_nQuad,1) + ...
                                     (B_2I_ig)'*(C)*(B_2I_ig)*M2N(idx_nQuad,1)); 
            end
    end
    % Body force
    b = f_b(xQuad(idx_nQuad,1),xQuad(idx_nQuad,2));
    FI_b = FI_b + (PSI')*b*(Weight(idx_nQuad,1));  
end % p_id end nodal integariton
sub_textprogressbar('done')    

%% Boundary Integration for Boundary Conditions
nQuad_BC = Quadrature.BC.nQuad_onBoundary;
xQuad_BC = Quadrature.BC.xQuad_onBoundary;
Weight_BC = Quadrature.BC.Weight_onBoundary;
Normal_BC = Quadrature.BC.Normal_onBoundary;
SHP_onBC = Quadrature.BC.SHP_BC;
SHPDX1_onBC = Quadrature.BC.SHPDX1_BC;
SHPDX2_onBC = Quadrature.BC.SHPDX2_BC;
sub_textprogressbar('Loop Over Quadrature Points on Boundaries, Progress: ')

for idx_nQuad = 1:nQuad_BC % for loop of quadrature points at boundary  
    sub_textprogressbar(100*idx_nQuad/nQuad_BC)      
        % normal at quadrature points 
        n1 = Normal_BC(idx_nQuad,1);
        n2 = Normal_BC(idx_nQuad,2);
        % PSI matrix at quadrature points 
        PSI(:,:) = 0;
        PSI(1,1:2:end) = SHP_onBC(idx_nQuad,:); 
        PSI(2,2:2:end) = SHP_onBC(idx_nQuad,:);
        switch Quadrature.BC.EBCtype{idx_nQuad}
            case {'NBC'}
            % surface traction
            t = f_t(xQuad_BC(idx_nQuad,1),xQuad_BC(idx_nQuad,2),n1,n2);
            FI_t = FI_t + PSI'*t*Weight_BC(idx_nQuad);
            case {'EBC'}
            % Surface normal Eta
            ETA = [n1 0 n2;
                   0 n2 n1];
            % Load B Matrix
            B(:,:) = 0; 
            B(1,1:2:end) = SHPDX1_onBC(idx_nQuad,:);
            B(2,2:2:end) = SHPDX2_onBC(idx_nQuad,:);
            B(3,1:2:end) = SHPDX2_onBC(idx_nQuad,:);
            B(3,2:2:end) = SHPDX1_onBC(idx_nQuad,:);
            % essential boundary condition g
            g = f_g(xQuad_BC(idx_nQuad,1),xQuad_BC(idx_nQuad,2));
            % Switch s
            S = f_S(xQuad_BC(idx_nQuad,1),xQuad_BC(idx_nQuad,2));
            % Nitche's term
            KIJ_g = KIJ_g + (PSI'*S*ETA*C*B)*Weight_BC(idx_nQuad);
            FI_g = FI_g + (B'*C*ETA'*S'*g)*Weight_BC(idx_nQuad);
            KIJ_beta = KIJ_beta + beta*(PSI'*S*PSI)*Weight_BC(idx_nQuad);
            FI_beta = FI_beta + beta*(PSI'*S*g)*Weight_BC(idx_nQuad);
            otherwise
        end % end if for boundary conditions
end % end for loop of quadrature points at boundary  
sub_textprogressbar('done')

%% Assembly
K_total = KIJ_c+KIJ_beta-(KIJ_g+KIJ_g');
F_total = FI_b+FI_t+FI_beta-FI_g;   
end
