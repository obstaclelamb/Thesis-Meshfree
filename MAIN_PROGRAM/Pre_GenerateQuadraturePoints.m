function [Quadrature] = Pre_GenerateQuadraturePoints(Quadrature,Discretization)
% This subroutine computes the quadrature points for domain integration and boundary integration for 
% nodal integration or Gauss integration.

% load nodal information, and quadrature information
xI = Discretization.xI; nP = Discretization.nP;
x1I = Discretization.xI(:,1); x2I = Discretization.xI(:,2);
x1I_Boundary = Discretization.xI_Boundary(:,1); x2I_Boundary = Discretization.xI_Boundary(:,2);
Index_BC = Discretization.Index_BC; % 1:length(x1I_Boundary);
nGaussPoints = Quadrature.nGaussPoints;
nGaussCells = Quadrature.nGaussCells;

switch Quadrature.Integration
    case {'SCNI','DNI'} % Nodal Integration method
        
        [VoronoiDiagram] = getVoronoiDiagram(x1I,x2I,Discretization);
        
        %% Nodal Integration
        % The number of evaluation poitns and their nodal position are
        nQuad = nP;  xQuad = xI;
        Weight = zeros(nQuad,1);
        for idx_quad = 1:nQuad
            Weight(idx_quad) = polyarea(VoronoiDiagram.VerticeCoordinates(VoronoiDiagram.VoronoiCell{idx_quad}(:),1),...
                                        VoronoiDiagram.VerticeCoordinates(VoronoiDiagram.VoronoiCell{idx_quad}(:),2));
        end
        
        Quadrature.VoronoiDiagram = VoronoiDiagram;
    otherwise % Gauss Integration
        %% Gauss Quadrature Points for Guass Integration
        % obtain the gauss points
        [nQuad,xQuad,Weight] = getBackgroundIntegrationCell(Discretization,xI,x1I_Boundary,x2I_Boundary,nGaussCells,nGaussPoints);
        % for Gauss integration, there is not voronoi information
%         VoronoiDiagram = struct([]);
        Quadrature.Table_XatBC = [];
end

% Save integration points for domain integration
Quadrature.Domain.xQuad = xQuad; 
Quadrature.Domain.nQuad = nQuad; 
Quadrature.Domain.Weight = Weight;

%% Generate quadrature points for boundary 
% Generate Gauss points for boundary Conditions
% Integrate over each line segment of the boundary points

switch Quadrature.Option_BCintegration % type of boundary condition
    case 'NODAL'
        Index_Vertices_onBoundary = VoronoiDiagram.Index_Vertices_onBoundary;
        xQuad_BC = []; Weight_BC = []; Normal_BC = [];
        VoronoiCell = VoronoiDiagram.VoronoiCell; % C
        XVoronoiVertices = VoronoiDiagram.VerticeCoordinates;  % V
%         V(C{ij},:) will give you the vertices of the ij'th data
        for idx_quad = 1:nQuad
        % Loop orver voronoi cell segment
            for k = 0:length(VoronoiCell{idx_quad})-1 % each edge for each cell
                % note that in voronoi cell, the vertices is label counter
                % clockwise, so the outnormal is always 90 degree
                if k == 0 % the last point of the k-th vertice of the cell
                    Vertex_onCell1 = VoronoiCell{idx_quad}(end);
                    x1Vertex_node1 = XVoronoiVertices(VoronoiCell{idx_quad}(end),1);
                    x2Vertex_node2 = XVoronoiVertices(VoronoiCell{idx_quad}(end),2);
                else
                    Vertex_onCell1 = VoronoiCell{idx_quad}(k);
                    x1Vertex_node1 = XVoronoiVertices(VoronoiCell{idx_quad}(k),1);
                    x2Vertex_node2 = XVoronoiVertices(VoronoiCell{idx_quad}(k),2);
                end
                Vertex_onCell2 = VoronoiCell{idx_quad}(k+1);
                x1Vertex_cell = XVoronoiVertices(VoronoiCell{idx_quad}(k+1),1);
                x2Vertex_cell = XVoronoiVertices(VoronoiCell{idx_quad}(k+1),2);
                
                % length and normal
                Lk_Cell = norm([x1Vertex_cell,x2Vertex_cell]-[x1Vertex_node1,x2Vertex_node2],2);
                x21_vertice = [x1Vertex_cell,x2Vertex_cell] - [x1Vertex_node1,x2Vertex_node2];
                VX21_NORMAL = x21_vertice*[cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
                Nk_Cell = VX21_NORMAL/(norm(VX21_NORMAL,2)+eps);
                if norm(Nk_Cell,2) < eps
                    Lk_Cell = eps;
                end
                
                MidpointXY_KthCellSegment = ([x1Vertex_node1 x2Vertex_node2] + [x1Vertex_cell x2Vertex_cell])/2;
                
                % Obtain the cell boundary length for the boundary
                % condition 
                if ismember(Vertex_onCell1,Index_Vertices_onBoundary) && ismember(Vertex_onCell2,Index_Vertices_onBoundary)
                xQuad_BC = [xQuad_BC; [MidpointXY_KthCellSegment]];
                Weight_BC = [Weight_BC; Lk_Cell];
                Normal_BC = [Normal_BC; [Nk_Cell]];
                end % end if for vertices are on boundary
            end
            
        end
    case 'GAUSS'
        N_Gauss_BC = Quadrature.nGaussPoints;
        xQuad_BC = []; Weight_BC = []; Normal_BC = [];
        for Meva_id = 1:length(Discretization.Index_BC)
            % two ends of each boudnary edge
            if Meva_id == length(Index_BC)  
            index_node1_BC = Index_BC(Meva_id);
            index_node2_BC = Index_BC(1);
            else
            index_node1_BC = Index_BC(Meva_id);
            index_node2_BC = Index_BC(Meva_id+1);
            end
            % generate Gauss points
            [x1Gauss_onBoundary, ~] = getGaussQuad(N_Gauss_BC, (xI([index_node1_BC],1)), (xI([index_node2_BC],1)));
            [x2Gauss_onBoundary, ~] = getGaussQuad(N_Gauss_BC, (xI([index_node1_BC],2)), (xI([index_node2_BC],2)));
            L_P21 = norm(xI(index_node2_BC,:)-xI(index_node1_BC,:),2);
            [~, WBCGAUSS] = getGaussQuad(N_Gauss_BC, 0, L_P21);
            % obtain normals and weights
            x21 = xI(index_node2_BC,:)-xI(index_node1_BC,:);
            x21_Normal = x21*[cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
            Normal_atBC = x21_Normal/(norm(x21_Normal,2)+eps);
 
            xQuad_BC = [xQuad_BC; [x1Gauss_onBoundary',x2Gauss_onBoundary']];
            Weight_BC = [Weight_BC; WBCGAUSS'];
            Normal_BC = [Normal_BC; [ones(length(WBCGAUSS),1)*Normal_atBC]];
            
        end % end NEVA for boundary integration points
        
        % if the interior hole exist, calculate those boundary as well
        if Discretization.DomainInteriorExist
            % find the Guass points within the domain due to outer boundary

            for idx_interior = 2:length(Discretization.xVertices_Inner)+1 % loop over all interior
                Index_BC_interior = Discretization.Index_BC_inner{idx_interior}(:);
%                 [idx_xInsideHole,~] = inpolygon(x1gauss_inside,x2gauss_inside,xI(Index_BC_interior,1),xI(Index_BC_interior,2));
                for Meva_id = 1:length(Index_BC_interior)
                    if Meva_id == length(Index_BC_interior)  
                    index_node1_BC = Index_BC_interior(Meva_id);
                    index_node2_BC = Index_BC_interior(1);
                    else
                    index_node1_BC = Index_BC_interior(Meva_id);
                    index_node2_BC = Index_BC_interior(Meva_id+1);
                    end
                    % generate Gauss points
                    [x1Gauss_onBoundary, ~] = getGaussQuad(N_Gauss_BC, (xI([index_node1_BC],1)), (xI([index_node2_BC],1)));
                    [x2Gauss_onBoundary, ~] = getGaussQuad(N_Gauss_BC, (xI([index_node1_BC],2)), (xI([index_node2_BC],2)));
                    L_P21 = norm(xI(index_node2_BC,:)-xI(index_node1_BC,:),2);
                    [~, WBCGAUSS] = getGaussQuad(N_Gauss_BC, 0, L_P21);
                    % obtain normals and weights
                    x21 = xI(index_node2_BC,:)-xI(index_node1_BC,:);
                    x21_Normal = x21*[cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
                    Normal_atBC = x21_Normal/(norm(x21_Normal,2)+eps);

                    xQuad_BC = [xQuad_BC; [x1Gauss_onBoundary',x2Gauss_onBoundary']];
                    Weight_BC = [Weight_BC; WBCGAUSS'];
                    Normal_BC = [Normal_BC; [ones(length(WBCGAUSS),1)*Normal_atBC]];
                end
            end
        end
        
end

NEVA_BC = length(Weight_BC);

% find the quadrature points which is EBC/BNC
EBCtype = cell(NEVA_BC,1);
EBCtype(:) = {'ZBC'}; % zero traction boundary condition
EBCtype(Discretization.CriteriaNBC(xQuad_BC(:,1),xQuad_BC(:,2))) = {'NBC'};
EBCtype(Discretization.CriteriaEBC(xQuad_BC(:,1),xQuad_BC(:,2))) = {'EBC'};

% save integration poitns on the boundary
Quadrature.BC.nQuad_onBoundary = NEVA_BC; 
Quadrature.BC.xQuad_onBoundary = xQuad_BC; 
Quadrature.BC.Weight_onBoundary = Weight_BC; 
Quadrature.BC.Normal_onBoundary = Normal_BC; % normal vector required for Nitches method
Quadrature.BC.EBCtype = EBCtype;

% save xI and # of node
Quadrature.xI = Discretization.xI; 
Quadrature.nP = Discretization.nP;


end

