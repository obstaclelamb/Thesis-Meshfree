function [SHPDX1_smoothed,SHPDX2_smoothed,Area_VoronoiCell] = getSmoothedDerivative(RK,xI,xV)
% This subroutine computes SCNI smoothed gradients for a given nodal representative domain
% INPUT PARAMETERS
%    RK:   structure contains all RK information
%    xI:   (RK nodes) vector NP by 2
%    xV:   (Cell vertices): Nodal representative doamin boundary vertices
%          coordintates

% OUTPUT PARAMETERS
%    SHPDX1_smoothed  - Smoothed gradient in x1
%    SHPDX2_smoothed  - Smoothed gradient in x2
%    Area_VoronoiCell - area of the nodal representative domain

nP = length(xI); % number of node  
nV = length(xV); % number of vertices  
% obtain the area of each voronoi cell
Area_VoronoiCell = polyarea(xV(:,1),xV(:,2))+eps;
% Loop orver voronoi cell segment
SHP_Smoothed_local = sparse(2,nP);
for k = 0:nV-1 % loop over each edge for each cell
    % find the two ends of the edge
    if k == 0 
        Vertex1 = xV(end,:);
    else
        Vertex1 = xV(k,:);
    end
    Vertex2 = xV(k+1,:);
    % Cell edge length and normal
    Lk_Cell = norm(Vertex2-Vertex1,2);
    xv21 = Vertex2 - Vertex1;
    xv21_Normal = xv21*[cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
    % calculate the unit normal, eps is to avoid singularity
    Nk_Cell = xv21_Normal/(norm(xv21_Normal,2)+eps);
    if norm(Lk_Cell,2) < eps || norm(Nk_Cell,2) < eps
        Lk_Cell = eps;
        Nk_Cell = [eps eps];
    end
    % Smoothing Points of each segment
    xtilde_CellEdge = (Vertex1 + Vertex2)/2;  
    % RK shape function evaluation at smoothing point
    [SHP_xtilde] = getRKShapeFunction(RK,xI,xtilde_CellEdge(1:2),[1,0,0]);
    SHP_Smoothed_local = SHP_Smoothed_local + (Nk_Cell'*SHP_xtilde)*Lk_Cell;
end % end of each voronoi cell boundary
% end of each voronoi cell boundary
SHPDX1_smoothed = (1/Area_VoronoiCell)*SHP_Smoothed_local(1,:);
SHPDX2_smoothed = (1/Area_VoronoiCell)*SHP_Smoothed_local(2,:);            
end