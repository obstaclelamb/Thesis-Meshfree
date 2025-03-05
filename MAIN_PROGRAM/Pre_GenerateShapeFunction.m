function [Quadrature] = Pre_GenerateShapeFunction(RK,Quadrature)
% This subroutine generates the RK shape functions and gradients at quadrature points.
%% COMPUTE RK SHAPE FUNCTION WITH DNI or SCNI
% load variables from quadrature points
nQuad = Quadrature.Domain.nQuad; % number of quadrature point
xQuad = Quadrature.Domain.xQuad; % quadrature point
nP = RK.nP; % number of node
xI = RK.xI; % node

% Load voronoi diagram
switch Quadrature.Integration
  case {'SCNI','DNI'}
    VoronoiDiagram = Quadrature.VoronoiDiagram;
end

SHP = sparse(nQuad,nP); % Shape Function
SHPDX1 = sparse(nQuad,nP); % Deravitive in x1
SHPDX2 = sparse(nQuad,nP); % Deravitive in x2
SHPDX1X1 = sparse(nQuad,nP); % Double Deravitive in x1
SHPDX2X2 = sparse(nQuad,nP); % Double Deravitive in x2
SHPDX1X2 = sparse(nQuad,nP); % Deravitive in x2 and x1

switch Quadrature.Integration
  case {'SCNI','DNI'}
    %% NODAL Integration
    % Evaluation of Shape Function
    xVoronoiVertices = VoronoiDiagram.VerticeCoordinates;
    VoroCELL = VoronoiDiagram.VoronoiCell;

    Area_VoronoiCell = zeros(nP,1); % Area of each cell
    Normal_VoronoiCell = cell(nP,1); % Normal of each cell
    SmoothingPoint_VoronoiCell = cell(nP,1); % Smoothing points of each cell
end
% Generate the shape function and it's deravitive at nodes for nodal integration

%% loop over integration points
sub_textprogressbar('Construction Shape Function, Progress:');
for idx_nQuad = 1:nQuad
  % progress bar
  sub_textprogressbar(100*idx_nQuad/nQuad)
  % determine the quadrature rule
  switch Quadrature.Integration
    case 'GAUSS'
      %% GAUSS Integration
      % Neighbor search
      [SHP(idx_nQuad,:),SHPDX1(idx_nQuad,:),SHPDX2(idx_nQuad,:)] = getRKShapeFunction(RK,xI,xQuad(idx_nQuad,:),[1,1,0]);

    case {'SCNI','DNI'}
      %% NODAL Integration
      % Evaluation of Shape Function
      xV = xVoronoiVertices(VoroCELL{idx_nQuad}(:),:);
      nV = length(xV(:,1));
      % obtain the area of each voronoi cell
      Area_VoronoiCell(idx_nQuad) = polyarea(xV(:,1),xV(:,2));
      switch Quadrature.Stabilization
        case {'N'} % N-type stabilization
          switch Quadrature.Integration
            case {'DNI'}
              [SHP(idx_nQuad,:),SHPDX1(idx_nQuad,:),SHPDX2(idx_nQuad,:),...
                SHPDX1X1(idx_nQuad,:),SHPDX2X2(idx_nQuad,:),SHPDX1X2(idx_nQuad,:)] = ...
                getRKShapeFunction(RK,xI,xQuad(idx_nQuad,:),[1,1,1],1);
            otherwise % SCNI
              [SHP(idx_nQuad,:),~,~,...
                SHPDX1X1(idx_nQuad,:),SHPDX2X2(idx_nQuad,:),SHPDX1X2(idx_nQuad,:)] = ...
                getRKShapeFunction(RK,xI,xQuad(idx_nQuad,:),[1,1,1],1);
          end
        otherwise % no stabilization, just evaluate shape function
          switch Quadrature.Integration
            case {'DNI'}
              [SHP(idx_nQuad,:),SHPDX1(idx_nQuad,:),SHPDX2(idx_nQuad,:)] = ...
                getRKShapeFunction(RK,xI,xQuad(idx_nQuad,:),[1,1,0]);
            case {'SCNI'} % SCNI
              [SHP(idx_nQuad,:)] = getRKShapeFunction(RK,xI,xQuad(idx_nQuad,:),[1,0,0]);
          end
      end

      SmoothedSHP_at_xtilde = zeros(2,nP);
      % USING MSCNI STABILIZATION
      switch Quadrature.Stabilization
        case {'M'} % M-type stabilization
          switch Quadrature.Integration
            case {'SCNI'} % MSCNI
              % Interior stabilization points
              Vertices1 = xV;
              Vertices2 = circshift(xV,[-1 0]);
              for s = 1 : nV % loop over the cell
                xV_SubCell = [Vertices1(s,:);
                  Vertices2(s,:);
                  xQuad(idx_nQuad,:);];

                [SHPDX1_Is{idx_nQuad}(s,:),SHPDX2_Is{idx_nQuad}(s,:),AREA_Is{idx_nQuad}(s,1)] = getSmoothedDerivative(RK,xI,xV_SubCell);
                if AREA_Is{idx_nQuad}(s,1) < eps
                  % for some cases, there will be additional ~0 area subcell, so we
                  % set the area and the corresponding shape function to be zero
                  SHPDX1_Is{idx_nQuad}(s,:) = zeros(1,nP);
                  SHPDX2_Is{idx_nQuad}(s,:) = zeros(1,nP);
                  AREA_Is{idx_nQuad}(s,1) = 0;
                end
              end
            case {'DNI'}
              % M-type Stabilization by cutting the smoothing zone into small pieces.
              xVerices_MSCNI = unique([xV; xQuad(idx_nQuad,:)],'rows');
              DT_MSCNI = delaunayTriangulation(xVerices_MSCNI(:,1),...
                xVerices_MSCNI(:,2)); % delauny for sub-triangle
              xNs_MSCNI = incenter(DT_MSCNI);
              [N_subcell,~] = size(xNs_MSCNI);
              for s = 1 : N_subcell % loop over the sub cell
                % vertices of the sub cell
                xV_SubCell = DT_MSCNI.Points(DT_MSCNI.ConnectivityList(s,:),:);
                [~,SHPDX1_Is{idx_nQuad}(s,:),SHPDX2_Is{idx_nQuad}(s,:)] = getRKShapeFunction(RK,xI,xNs_MSCNI(s,:),[1,1,0]);
                AREA_Is{idx_nQuad}(s,1) = polyarea(xV_SubCell(:,1),xV_SubCell(:,2));
              end
          end

        case {'N'}
          % NSNI stabilization
          Vertices1 = xV - ones(nV,1)*mean(xV);
          Vertices2 = circshift(xV - ones(nV,1)*mean(xV),[-1 0]);
          % using formula to calculate innertia for arbitrary polygon
          M1N(idx_nQuad) = (1/12)*abs(sum((Vertices1(:,2).^2+Vertices1(:,2).*Vertices2(:,2)+Vertices2(:,2).^2)...
            .*(Vertices1(:,1).*Vertices2(:,2) - Vertices2(:,1).*Vertices1(:,2))));
          M2N(idx_nQuad) = (1/12)*abs(sum((Vertices1(:,1).^2+Vertices1(:,1).*Vertices2(:,1)+Vertices2(:,1).^2)...
            .*(Vertices1(:,1).*Vertices2(:,2) - Vertices2(:,1).*Vertices1(:,2))));
          dx1 = sqrt((xQuad(idx_nQuad,1)-mean(xV(:,1))).^2);
          dx2 = sqrt((xQuad(idx_nQuad,2)-mean(xV(:,2))).^2);
          % obtain intertia for NSNI
          M1N(idx_nQuad) = M1N(idx_nQuad) + Area_VoronoiCell(idx_nQuad)*dx1^2;
          M2N(idx_nQuad) = M2N(idx_nQuad) + Area_VoronoiCell(idx_nQuad)*dx2^2;
      end

      % Loop orver voronoi cell segment
      for k = 0:nV-1 % each edge for each cell
        % note that in voronoi cell, the vertices is label counter
        % clockwise, so the outnormal is always 90 degree
        if k == 0 % the last point of the k-th vertice of the cell
          Vertex1 = xV(end,:);
        else
          Vertex1 = xV(k,:);
        end
        Vertex2 = xV(k+1,:);

        % length and normal
        Lk_Cell = norm(Vertex2-Vertex1,2);
        xv21 = Vertex2 - Vertex1;
        xv21_Normal = xv21*[cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
        % calculate the unit normal, eps is to avoid singularity
        Nk_Cell = xv21_Normal/(norm(xv21_Normal,2)+eps);
        if norm(Nk_Cell,2) < eps
          Lk_Cell = eps;
        end
        % Smoothing Points of each segment
        xtilde_CellEdge = (Vertex1 + Vertex2)/2;
        Normal_VoronoiCell{idx_nQuad} = [Normal_VoronoiCell{idx_nQuad},Nk_Cell'];
        SmoothingPoint_VoronoiCell{idx_nQuad} = [SmoothingPoint_VoronoiCell{idx_nQuad},xtilde_CellEdge'];

        %% Here we need to formulate the smoothed deravitive for SCNI
        switch Quadrature.Integration
          case {'SCNI'}
            % RK shape function
            [SHP_Smoothed] = getRKShapeFunction(RK,xI,xtilde_CellEdge(1:2),[1,0,0]);
            SmoothedSHP_at_xtilde = SmoothedSHP_at_xtilde + (Nk_Cell'*SHP_Smoothed)*Lk_Cell;
            %                 SHP_LOCAL = SHP_LOCAL + norm(MidpointXY_CELL-XI(p_id),2)*SHPD_XI;
            if sum(isnan(Nk_Cell))>0
              stop
            end

        end
      end % end of each voronoi cell boundary
      %% Here we need to formulate the smoothed deravitive for SCNI
      switch Quadrature.Integration
        case {'SCNI'} % SCNI
          SHPDX1(idx_nQuad,:) = (1/Area_VoronoiCell(idx_nQuad))*SmoothedSHP_at_xtilde(1,:);
          SHPDX2(idx_nQuad,:) = (1/Area_VoronoiCell(idx_nQuad))*SmoothedSHP_at_xtilde(2,:);
      end
  end % end switch integration method
end % end p_id loop for constructing shape function at RK nodes
% delete(WB);
sub_textprogressbar('done')

%% save information for stabilization of MSCNI/NSNI
switch Quadrature.Integration
  case {'SCNI','DNI'}
    switch Quadrature.Stabilization
      case {'M'}
        Quadrature.MtypeStablization.SHPDX1_Is = SHPDX1_Is;
        Quadrature.MtypeStablization.SHPDX2_Is = SHPDX2_Is;
        Quadrature.MtypeStablization.AREA_Is = AREA_Is;

      case {'N'}
        % The smoothed second order deravitive does not work
        Quadrature.NtypeStablization.SHPDX1X1 = SHPDX1X1;
        Quadrature.NtypeStablization.SHPDX2X2 = SHPDX2X2;
        Quadrature.NtypeStablization.SHPDX1X2 = SHPDX1X2;
        Quadrature.NtypeStablization.SHPDX2X1 = SHPDX1X2;
        Quadrature.NtypeStablization.M = zeros(nP,2);
        Quadrature.NtypeStablization.M(:,1) = M1N;
        Quadrature.NtypeStablization.M(:,2) = M2N;
    end
end

%% evaluating shape function at boundary points
nQuad_BC = Quadrature.BC.nQuad_onBoundary;
xQuad_BC = Quadrature.BC.xQuad_onBoundary;
SHP_BC = zeros(nQuad_BC,nP);
SHPDX1_BC = zeros(nQuad_BC,nP);
SHPDX2_BC = zeros(nQuad_BC,nP);
sub_textprogressbar('Construction Shape Function at Boundary, Progress:');
for idx_Neva_BC = 1:nQuad_BC
  sub_textprogressbar(100*idx_Neva_BC/nQuad_BC)
  [SHP_BC(idx_Neva_BC,:),SHPDX1_BC(idx_Neva_BC,:),SHPDX2_BC(idx_Neva_BC,:)] = getRKShapeFunction(RK,xI,xQuad_BC(idx_Neva_BC,:),[1,1,0]);
end
sub_textprogressbar('done');


%% Test Reproducing Condition of the Shape Function and Its Deravitive
[PASS_ORNOT] = sub_TestCompleteness(RK.Order,xI,xQuad,SHP,SHPDX1,SHPDX2,'OFF');

if PASS_ORNOT
  disp('RK Shape function pass reproducing condition')
end

%% Store the shape function in RK structure as well in favor of future use
Quadrature.SHP = SHP;
Quadrature.SHPDX1 = SHPDX1;
Quadrature.SHPDX2 = SHPDX2;

Quadrature.BC.SHP_BC = SHP_BC;
Quadrature.BC.SHPDX1_BC = SHPDX1_BC;
Quadrature.BC.SHPDX2_BC = SHPDX2_BC;

end

