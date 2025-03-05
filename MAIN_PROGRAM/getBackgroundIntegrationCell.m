function [nQuad,xQuad,Weight] = getBackgroundIntegrationCell(Discretization,xI,x1I_Boundary,x2I_Boundary,nGaussCells,nGaussPoints)
% This subroutine constructs the rectangular background Gauss integration cells if Gauss integration is adopted.
% Generate Background integration cell and points
% check the shorter side of the problem domain, is in x or y?
x1GaussCell_domainsize = norm(max(xI(:,1))-min(xI(:,1))); 
x2GaussCell_domainsize = norm(max(xI(:,2))-min(xI(:,2)));
nGauss_in_x2_direction = fix(nGaussCells*x2GaussCell_domainsize/min([x1GaussCell_domainsize,x2GaussCell_domainsize]));
nGauss_in_x1_direction = fix(nGaussCells*x1GaussCell_domainsize/min([x1GaussCell_domainsize,x2GaussCell_domainsize]));
nGauss_in_x2_direction = max([nGauss_in_x2_direction,nGaussCells]);
nGauss_in_x1_direction = max([nGauss_in_x1_direction,nGaussCells]);
x1_backgroundcell = linspace(min(xI(:,1)), max(xI(:,1)),nGauss_in_x1_direction+1);
x2_backgroundcell = linspace(min(xI(:,2)), max(xI(:,2)),nGauss_in_x2_direction+1);
x1gauss = []; x2guass = []; weight_gauss = [];
% for given cell, generate the Gauss point in x and y direction
for iX = 1:length(x1_backgroundcell)-1
    for iY = 1:length(x2_backgroundcell)-1
        [X1GAUSS, WX1GAUSS] = getGaussQuad(nGaussPoints, x1_backgroundcell(iX), x1_backgroundcell(iX+1));
        [X2GAUSS, WX2GAUSS] = getGaussQuad(nGaussPoints, x2_backgroundcell(iY), x2_backgroundcell(iY+1));
        nGaussPerCell = nGaussPoints*nGaussPoints;
        [x1G_unsort,x2G_unsort] = meshgrid(X1GAUSS,X2GAUSS);
        [W_x1G,W_X2G] = meshgrid(WX1GAUSS,WX2GAUSS);
        x1gauss = [x1gauss;reshape(x1G_unsort,[nGaussPerCell 1])];
        x2guass = [x2guass;reshape(x2G_unsort,[nGaussPerCell 1])];
        weight_gauss = [weight_gauss;reshape(W_x1G.*W_X2G,[nGaussPerCell 1])]; 
    end
end


if Discretization.DomainInteriorExist
    % find the Guass points within the domain due to outer boundary
    [idx_xInside,~] = inpolygon(x1gauss,x2guass,x1I_Boundary,x2I_Boundary);
    x1gauss_inside = x1gauss(idx_xInside); x2gauss_inside = x2guass(idx_xInside);
    weight_gauss_inside = weight_gauss(idx_xInside); nQuad = length(weight_gauss_inside);
    
    for idx_interior = 2:length(Discretization.xVertices_Inner)+1 % loop over all interior
        Index_BC_interior = Discretization.Index_BC_inner{idx_interior}(:);
        [idx_xInsideHole,~] = inpolygon(x1gauss_inside,x2gauss_inside,xI(Index_BC_interior,1),xI(Index_BC_interior,2));
        x1gauss_inside = x1gauss_inside(~idx_xInsideHole); x2gauss_inside = x2gauss_inside(~idx_xInsideHole);
        weight_gauss_inside = weight_gauss_inside(~idx_xInsideHole); 
        nQuad = length(weight_gauss_inside);
    end
else
    % find the Guass points within the domain
    [idx_xInside,~] = inpolygon(x1gauss,x2guass,x1I_Boundary,x2I_Boundary);
    x1gauss_inside = x1gauss(idx_xInside); x2gauss_inside = x2guass(idx_xInside);
    weight_gauss_inside = weight_gauss(idx_xInside); nQuad = length(weight_gauss_inside);
end

% display
disp(['Gauss Integration: the background quadratures points: ',num2str(nQuad),' points'])

xQuad = [x1gauss_inside,x2gauss_inside];
Weight = weight_gauss_inside;
% for Gauss integration, there is not voronoi information
% VoronoiDiagram = struct([]);


end

