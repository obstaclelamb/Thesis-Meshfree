function [VORONOID] = getVoronoiDiagram(x1I,x2I,Discretization)
    %% Genearating Voronoi cell with boundary segmentation
    % This subroutine calls sub_VoronoiLimit.m to generate Voronoi cells and 
    % then re-arranges the Voronoi cell IDs and coordinates.
    % Note: There are several posibility that you cannot obtain the voronoi
    % cell with the prescribed boundary: 
    % 1. The domain geometry has too strong convex
    % 2. Boundary points are too dense and interior points are too sparse
    % Some geometry will cause problem inside the delauny diagram in the sub_VoronoiLimit.
    % The good way to avoid it is to avoid bad discretization 
    
    % inner-domain
    xP = [x1I,x2I];
    if Discretization.DomainInteriorExist % if there is any interior hole exist
    InteriorCoordinatesCell = Discretization.xI_InteriorBoundary;
    [VerticeCoordinates,VoronoiCELL_OLD,xI_new,Vertices_onBoundary]  = sub_VoronoiLimit(x1I,x2I,'bs_ext',[Discretization.xI_Boundary],...
                                                                                             'bs_int',InteriorCoordinatesCell,'figure','off');
    else % no interior hole
        [VerticeCoordinates,VoronoiCELL_OLD,xI_new,Vertices_onBoundary]  = sub_VoronoiLimit(x1I,x2I,'bs_ext',[Discretization.xI_Boundary],'figure','off');
    
    end
    % Since the package will re-order the x1 and x2 coordinate, I have to re-order it
    % again to make sure the index is identical to the origianl one
    [~,index_XP] = sortrows(xP); [~,index_XI] = sortrows(xI_new);
    NODAL_INDEX(index_XP) = index_XI;
    VoronoiCell(:) = {VoronoiCELL_OLD{NODAL_INDEX(:)}};
 
    % Also, we need to find out the vertices on the boundary
    % For very rare situation, it cannot find all vertices, may be due to the
    % error in the coordinate are larger than precision of kn search
    [Index_Vertices_onBoundary,~] = knnsearch(VerticeCoordinates,Vertices_onBoundary,'K',1);
    VORONOID.Index_Vertices_onBoundary = Index_Vertices_onBoundary;

    VORONOID.VerticeCoordinates = VerticeCoordinates;
    VORONOID.VoronoiCell = VoronoiCell;
    VORONOID.Vertices_onBoundary = Vertices_onBoundary;
end

