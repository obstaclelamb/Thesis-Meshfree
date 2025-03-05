function [Discretization,Model] = Pre_GenerateDiscretization(varargin)
% This subroutine generates point-based domain discretization based on 
% domain vertices coordinates and the chosen discretization method
%
% Input: Model structure with the following information
%        xVertices
% Output: Discretization: Discretization structure define the RK nodal position

% if only Model is imported, then use default mesh algorithm by Matlab
Model = varargin{1};

% outer boundary
x1_vertices = Model.xVertices(:,1);
x2_vertices = Model.xVertices(:,2);

% for the version of short course, the interior hole is assumed cannot be exist so 
% we always set it to zero, if one want to model the problem, one can turn
% this function = 1 
Model.DomainInteriorExist = 0;
% % and set this hidden card in the input files to define the geometry of the
% % interior domain
% Model.DomainInteriorExist = 1; % if inner hole exist
% % create the inner domain for Voronoi tesselation, any polygon works for
% Model.Interior.NumberOfInterior = 2;
% % Create the vertices for the interior
% % format: cell{i} means the i-th interior hole
% % inside the cell{i} is the vertices coordinate the that hole
% % e,g, xVertices_Inner{1} = [[1.2 4.3 4 1]', [1.9 2 4.5 4]']; 
% xVertices_Inner = cell(Model.Interior.NumberOfInterior,1);
% xVertices_Inner{1} = [[1.2 4.3 4 1]', [1.9 2 4.5 4]'];
% xVertices_Inner{2} = [[5 7 8 7 5]', [1 0.5 4.5 5.5 4]'];
% % xVertices_Inner{3} = [[8.5 9.5 9.5 8.5]', [0.3 0.8 0.6 0.3]'];
% Model.xVertices_Inner = xVertices_Inner;
% % p.s only discretization method 'A' works for the this hidden card


% this new function is used for any new interior hole
if Model.DomainInteriorExist % if there is any interior hole exist
    DomainArea = polyarea(x1_vertices,x2_vertices); % outer boundary domain
    for idx_interior = 1:length(Model.xVertices_Inner) % loop over all interior
    x1_vertices_Inner = Model.xVertices_Inner{idx_interior}(:,1);
    x2_vertices_Inner = Model.xVertices_Inner{idx_interior}(:,2);
    DomainArea = DomainArea - polyarea(x1_vertices_Inner,x2_vertices_Inner);
    end
else
    DomainArea = polyarea(x1_vertices,x2_vertices);
end


switch Model.Discretization.Method 
    case {'A'}
        % create PDE model object for MATLAB default mesh generator
        % for more information, please refer to 
        % https://www.mathworks.com/help/pde/ug/pde.pdemodel.geometryfromedges.html?s_tid=doc_ta
        model = createpde;
        % read in the domain vertices from Model of getInput
        x1_vertices = Model.xVertices(:,1);
        x2_vertices = Model.xVertices(:,2);
        % Create polygon based on edge of the vertices
        R1 = [3,length(x1_vertices),x1_vertices',x2_vertices']';
        
        %% if there is a (or more hole), using this subroutine to generate
        % mesh accordingly from the description in the begining of this
        % subroutine 
        if Model.DomainInteriorExist % if there is any interior hole exist
            Discretization.DomainInteriorExist = 1;
            % if 50 exist the max vertices of the domain, just + it.
            R_AllBoundary = zeros(10,1+length(Model.xVertices_Inner));
            R_AllBoundary(1:length(R1),1) = R1;
            
            % it is a bad idea to use eval, but I have no choice in order
            % to automatically generate the input geometry...
            eval(['R' num2str(1),'=R_AllBoundary(:,1);']);
            % sf is a set formula created to define if there is any subtraction 
            % between geometry, eg sf = R1-C1 where C1 may be a circle  
            sf = ['R1']; % initialize the really outer domain
            ns = ['R1']; % initialize the really outer domain
            % loop over each hole domain
            for idx_interior = 1:length(Model.xVertices_Inner) % loop over all interior
            x1_vertices_Inner = Model.xVertices_Inner{idx_interior}(:,1);
            x2_vertices_Inner = Model.xVertices_Inner{idx_interior}(:,2);
            R_AllBoundary(1:2+2*length(x2_vertices_Inner),idx_interior+1) = [3,length(x1_vertices_Inner),x1_vertices_Inner',x2_vertices_Inner']';
            % it is a bad idea to use eval, but I have no choice in order
            % to automatically generate the input geometry...
            eval(['R' num2str(idx_interior+1),'=R_AllBoundary(:,idx_interior+1);']);
            % sf is a set formula created to define if there is any subtraction 
            % between geometry, eg sf = R1-C1 where C1 may be a circle  
            sf = [sf,'-R',num2str(idx_interior+1)];
            % create geometry
            ns = [ns; ['R',num2str(idx_interior+1)];];
            end

            % gm is the geometry based on vertices, 
            % sf is a set formula created to define if there is any subtraction 
            % between geometry, eg sf = R1-C1 where C1 may be a circle 
            
            gm = R_AllBoundary; 
            ns = ns';
            % Decompose constructive geometry into minimal regions 
            g = decsg(gm,sf,ns);
            % create geometry for Model object 
            geoModel = geometryFromEdges(model,g);
            % generate FE mesh for Model by MATLAB
            %         
            if Model.Discretization.Hmax <= 0
                FEmesh = generateMesh(model,'GeometricOrder','linear');
            else
                FEmesh = generateMesh(model,'Hmax',Model.Discretization.Hmax,'GeometricOrder','linear');
            end
            % Use FE mesh nodes as RK nodes.
            RKnodes = FEmesh.Nodes';  

            x1I = RKnodes(:,1);
            x2I = RKnodes(:,2);

            % find out the nodes in and on the domain
            [in_idx,on_idx] = inpolygon(x1I,x2I,x1_vertices,x2_vertices);

            % find notes that is not on the interior boundaries
            total_on_idx_inner = 0*[1:length(x1I)]';
%             total_on_idx_inner = [];
            xI_InteriorBoundary = cell(length(Model.xVertices_Inner),1);
            for idx_interior = 1:length(Model.xVertices_Inner) % loop over all interior
            x1_vertices_Inner = Model.xVertices_Inner{idx_interior}(:,1);
            x2_vertices_Inner = Model.xVertices_Inner{idx_interior}(:,2);
            [~,on_idx_inner] = inpolygon(x1I,x2I,x1_vertices_Inner,x2_vertices_Inner); % inner-domain
            
            % also, construct the boundary nodes for the interior hoel
            xV_Inner = [x1_vertices_Inner; x1_vertices_Inner(1)];
            yV_Inner = [x2_vertices_Inner; x2_vertices_Inner(1)];
            polygonedge_length = sum(sqrt(gradient(xV_Inner).^2 + gradient(yV_Inner).^2));
            nP_on_edges = length(find(on_idx_inner>0));
            
            % save the interior domain in a cell form
            [x1I_Boundary_inner,x2I_Boundary_inner] = interpm(xV_Inner,yV_Inner,polygonedge_length/nP_on_edges);
            x1I_Boundary_inner(end) = []; x2I_Boundary_inner(end) = [];
            xI_InteriorBoundary{idx_interior} = [x1I_Boundary_inner,x2I_Boundary_inner];

%             total_on_idx_inner = setdiff(total_on_idx_inner,find(on_idx_inner));
            total_on_idx_inner = total_on_idx_inner + on_idx_inner;
            end
%             idx_realinter = [find(in_idx&~on_idx);find(~on_idx_inner)];
%             idx_realinter

            total_on_idx_inner(total_on_idx_inner>1) = 1;
            % nodes in the interior
            x1I_Interior = x1I(in_idx&~on_idx&~total_on_idx_inner); 
            x2I_Interior = x2I(in_idx&~on_idx&~total_on_idx_inner);
            
            % really outer boundary
            xV = [x1_vertices; x1_vertices(1)];
            yV = [x2_vertices; x2_vertices(1)];
            polygonedge_length = sum(sqrt(gradient(xV).^2 + gradient(yV).^2));
            nP_on_edges = length(find(on_idx>0));

            % outer domain boundary
            [x1I_Boundary,x2I_Boundary] = interpm([x1_vertices; x1_vertices(1)],[x2_vertices; x2_vertices(1)],polygonedge_length/nP_on_edges);
            x1I_Boundary(end) = []; x2I_Boundary(end) = [];

        %% for domain with only outer boundary
        else 
            % no interior domain
            Discretization.DomainInteriorExist = 0;
            
            % gm is the geometry based on vertices, 
            % sf is a set formula created to define if there is any subtraction 
            % between geometry, eg sf = R1-C1 where C1 may be a circle  
            gm = [R1]; sf = 'R1';
            % create geometry
            ns = char('R1'); ns = ns';
            % Decompose constructive geometry into minimal regions 
            g = decsg(gm,sf,ns);
            % create geometry for Model object 
            geoModel = geometryFromEdges(model,g);
            % generate FE mesh for Model by MATLAB
            %         
            if Model.Discretization.Hmax <= 0
                FEmesh = generateMesh(model,'GeometricOrder','linear');
            else
                FEmesh = generateMesh(model,'Hmax',Model.Discretization.Hmax,'GeometricOrder','linear');
            end
            % Use FE mesh nodes as RK nodes.
            RKnodes = FEmesh.Nodes';  

            x1I = RKnodes(:,1);
            x2I = RKnodes(:,2);

            % find out the nodes in and on the domain
            [in_idx,on_idx] = inpolygon(x1I,x2I,x1_vertices,x2_vertices);
            x1I_Interior = x1I(in_idx&~on_idx); x2I_Interior = x2I(in_idx&~on_idx);

            xV = [x1_vertices; x1_vertices(1)];
            yV = [x2_vertices; x2_vertices(1)];
            polygonedge_length = sum(sqrt(gradient(xV).^2 + gradient(yV).^2));
            nP_on_edges = length(find(on_idx>0));

            [x1I_Boundary,x2I_Boundary] = interpm([x1_vertices; x1_vertices(1)],[x2_vertices; x2_vertices(1)],polygonedge_length/nP_on_edges);
            x1I_Boundary(end) = []; x2I_Boundary(end) = [];
        end
  
    case {'B'}
        %% if one want to use random mesh 
        % manually define the discretization by oneself
        % This is one example that we generate some background nodes
        % and only used the nodes within the domain
        nx1 = Model.Discretization.nx1 ;     % Number of node in x
        nx2 = Model.Discretization.nx2 ;     % Number of node in y
        
        np_backgroundmesh = (nx1)*(nx2);

        xc= linspace(min(x1_vertices),max(x1_vertices),nx1);
        dx = max(diff(xc));
        yc= linspace(min(x2_vertices),max(x2_vertices),nx2);
        dy = max(diff(xc));
        
        [x,y] = meshgrid(xc,yc);
        x1I = reshape(x,[(nx1)*(nx2) 1]);
        x2I = reshape(y,[(nx1)*(nx2) 1]);
        % random distribution of the nodes

        % find out the nodes in and on the polygon
        [in_idx,on_idx] = inpolygon(x1I,x2I,x1_vertices,x2_vertices);
        x1I_Interior = x1I(in_idx&~on_idx); x2I_Interior = x2I(in_idx&~on_idx);
        
        % perturb the interior by randomness
        x1I_Interior  = x1I_Interior + Model.Discretization.Randomness*dx*(rand(size(x1I_Interior))-0.5);
        x2I_Interior  = x2I_Interior + Model.Discretization.Randomness*dy*(rand(size(x2I_Interior))-0.5);
      
        
        % The nodes on the boundary are also computed manually
        x1I_Boundary = flipud([x(1:end,1);x(end,2:end)';x(end-1:-1:1,end);x(1,end-1:-1:2)']);
        x2I_Boundary = flipud([y(1:end,1);y(end,2:end)';y(end-1:-1:1,end);y(1,end-1:-1:2)']);
        
        % perturb the boundary by randomness
        target1 = find(x1I_Boundary<min(x1_vertices)+1E-7 & x2I_Boundary>min(x2_vertices)+1E-7 & x2I_Boundary<max(x2_vertices)-1E-7);
        target2 = find(x1I_Boundary>max(x1_vertices)-1E-7 & x2I_Boundary>min(x2_vertices)+1E-7 & x2I_Boundary<max(x2_vertices)-1E-7);
        target3 = find(x2I_Boundary<min(x2_vertices)+1E-7 & x1I_Boundary>min(x1_vertices)+1E-7 & x1I_Boundary<max(x1_vertices)-1E-7);
        target4 = find(x2I_Boundary>max(x2_vertices)-1E-7 & x1I_Boundary>min(x1_vertices)+1E-7 & x1I_Boundary<max(x1_vertices)-1E-7);
%         
        x2I_Boundary(target1) = x2I_Boundary(target1) + Model.Discretization.Randomness*(dy)*(rand(size(x2I_Boundary(target1)))-0.5);
        x2I_Boundary(target2) = x2I_Boundary(target2) + Model.Discretization.Randomness*(dy)*(rand(size(x2I_Boundary(target2)))-0.5);
        x1I_Boundary(target3) = x1I_Boundary(target3) + Model.Discretization.Randomness*(dx)*(rand(size(x2I_Boundary(target3)))-0.5);
        x1I_Boundary(target4) = x1I_Boundary(target4) + Model.Discretization.Randomness*(dx)*(rand(size(x2I_Boundary(target4)))-0.5);   
        
    if Model.DomainInteriorExist
        error('The domain contains a hole, where only discretization method A works, please turn to A');
    end
    case {'C'}
        %% if one want to use distorted mesh for the plate problem
        % based on the original Fortran program by Alek Shestakov; that program 
        % is found in: Nuclear Science & Engineering, Vol 105, pp.88-104 (1990),
        % "Test Problems in Radiative Transfer Calculations", by A. Shestakov, 
        nc = Model.Discretization.nc; 
        distortion = Model.Discretization.Distortion;
        [z0,r0] = sub_Shestakov(nc,distortion);

        % map the distorted mesh domain to the plate with a hole domain
        ztheta1 = (pi/4)*z0;
        rtheta1 = (3*r0 + 1) + 4*(sqrt((tan(ztheta1)).^2 + 1^2) - (1)).*r0;
        [x_mesh1,y_mesh1] = pol2cart(ztheta1,rtheta1);

        [z0,r0] = sub_Shestakov(nc,distortion);
        ztheta1 = (-pi/4)*z0 ;
        rtheta1 = (3*r0 + 1) + 4*(sqrt((tan(ztheta1)).^2 + 1^2) - (1)).*r0;
        rtheta1 = (rtheta1); 
        ztheta1 = (ztheta1)+pi/2; 
        
        
        [x_mesh2,y_mesh2] = pol2cart(ztheta1,rtheta1);

        xI_mesh = [flipud(x_mesh2(2:end,:)); x_mesh1]; yI_mesh = [flipud(y_mesh2(2:end,:)); y_mesh1];
        % figure(1), plot(x_mesh1,y_mesh1,'k-',x_mesh1',y_mesh1','k-','LineWidth',1.5); hold on;
        % figure(1), plot(x_mesh2,y_mesh2,'r-',x_mesh2',y_mesh2','r-','LineWidth',1.5); hold on;
        % figure, plot(xI_mesh,yI_mesh,'r-',xI_mesh',yI_mesh','r-','LineWidth',1.5); hold on;
        [theta_mesh,r_mesh] = cart2pol(xI_mesh,yI_mesh);
        % figure, surf(xI_mesh,yI_mesh);

        theta_edge = [theta_mesh(end,1:end-1) flip(theta_mesh(2:end,end)') flip(theta_mesh(1,2:end)) (theta_mesh(1:end-1,1)')]';
        r_edge = [r_mesh(end,1:end-1) flip(r_mesh(2:end,end)') flip(r_mesh(1,2:end)) (r_mesh(1:end-1,1)')]';
        theta_edge(end) = []; r_edge(end) = [];

        theta_INTERIOR = reshape(theta_mesh(2:end-1,2:end-1),[numel(theta_mesh(2:end-1,2:end-1)) 1]);
        r_INTERIOR = reshape(r_mesh(2:end-1,2:end-1),[numel(r_mesh(2:end-1,2:end-1)) 1]);

        [x1I_Boundary,x2I_Boundary] = pol2cart(theta_edge,r_edge);
        [x1I_Interior,x2I_Interior] = pol2cart(theta_INTERIOR,r_INTERIOR);
    if Model.DomainInteriorExist
        error('The domain contains a hole, where only discretization method A works, please turn to A');
    end
    case {'D'}
        %% if one want to use CAD model for the discretization
        
        [xI_fromCAD,BC_index] = sub_ReadNeutralInputFiles(Model.Discretization.InputFileName);
        
        % once you read in the coordinates, you need to make the boundary
        % nodal index follow the correct order (CCW + close loop). The
        % "unstructured" boundary index will lead to wrong 

        % The specifies shrink factor s using 
        % any of the previous syntaxes. s is a scalar between 0 and 1. 
        % Setting s to 0 gives the convex hull, and 
        % setting s to 1 gives a compact boundary that envelops the points. 
        % The default shrink factor is 0.5.
        
        s = 0.5; % shrink factor
        % ===================================================
        % this one sometimes give reated index, have not resolve this issue
        % yet. Hence, a switch here is employ such that if repeated index
        % appear, then use the stable one
%         wew = size(xI_fromCAD,2) - 1;
%         [~,k]= find_2or3d_boundary(xI_fromCAD,wew);
% %         stop
%         k = [];
%         while s < 1.001
%         [k] = boundary([xI_fromCAD(BC_index,1:2)],s);
%         idx_bdy = BC_index(k(1:end-1));
%         xI_fromCAD_bdy_unique = unique(xI_fromCAD(idx_bdy,:),'rows','stable');
% %         k = unique(k);
%         if (size(xI_fromCAD_bdy_unique) == size(xI_fromCAD(BC_index,:)))
%             break
%         end
%         s = s + 0.01;
%         end
        % ===================================================
        
        while s < 1.0
        [k] = boundary([xI_fromCAD(:,1:2)],s);
        idx_bdy = k(1:end-1);
        xI_fromCAD_bdy_unique = unique(xI_fromCAD(idx_bdy,:),'rows','stable');
%         k = unique(k);
        if (size(xI_fromCAD_bdy_unique) == size(xI_fromCAD(BC_index,:)))
        if (sort(xI_fromCAD_bdy_unique) == sort(xI_fromCAD(BC_index,:)))
%         if (size(xI_fromCAD_bdy_unique) == size(xI_fromCAD(BC_index,:)))
            break
        end
        end
        s = s + 0.01;
        end
        
        % if s > 1, that means the above method does not work
        if s > 1
%         % ===================================================
%         % this one is stable, but cannot hand convex geometry
        [k] = boundary(xI_fromCAD(:,1), xI_fromCAD(:,2),0.92);
        idx_bdy = k(1:end-1);
%         % ===================================================
        end
        
        x1I_Boundary = xI_fromCAD(idx_bdy,1); x2I_Boundary = xI_fromCAD(idx_bdy,2);
        x1I_Interior = xI_fromCAD(setdiff(1:end,idx_bdy),1); 
        x2I_Interior = xI_fromCAD(setdiff(1:end,idx_bdy),2);
        
        % this line is plotting for debug
%         plot(xI_fromCAD(:,1),xI_fromCAD(:,2),'ko',x1I_Interior, x2I_Interior,'b*',x1I_Boundary, x2I_Boundary,'gV-'); hold on


    if Model.DomainInteriorExist
        error('The domain contains a hole, where only discretization method A works, please turn to A');
    end
end



% Define the boundaries to be counter clockwise !
% This step is important for generating Voronoid diagram
[x1I_Boundary, x2I_Boundary] = poly2ccw(x1I_Boundary, x2I_Boundary);

if ispolycw(x1I_Boundary,x2I_Boundary)
    error('The boundary nodes are not following counter-clockwise direction, please re-check them');
    plot(x1I_Boundary,x2I_Boundary,'ro-'); hold on;
end

% Locate the information 
Discretization.xVertices = [x1_vertices, x2_vertices];
Discretization.xI_Boundary = [x1I_Boundary, x2I_Boundary];
Discretization.xI_Interior = [x1I_Interior, x2I_Interior];
Discretization.xI = [Discretization.xI_Boundary; Discretization.xI_Interior;];
Discretization.Area_Poly = DomainArea;

% Number of Boundary Points

% Assign the boundary point's index, or can be generated by input file
% Note: the index for each points should follow counter-clockwise direction
% for the polygon geometery

%% inner-domain
% for domain with inner hole
if Model.DomainInteriorExist
    % for this version, the domain interior is not permitted
    Discretization.DomainInteriorExist = 1;
    Discretization.xI_InteriorBoundary = xI_InteriorBoundary;

    N_BC = length(x1I_Boundary);
    N_BC_Interior = sum(cellfun(@length,xI_InteriorBoundary));
    % x1I = [x1I_Boundary; x1I_Boundary_inner; x1I_Interior;]; 
    % x2I = [x2I_Boundary; x2I_Boundary_inner; x2I_Interior;];
    x1I = [x1I_Boundary;]; 
    x2I = [x2I_Boundary;];
    Index_BC_inner = cell(length((Model.xVertices_Inner))+1,1);
    Index_BC_inner{1} = 1:N_BC;
    for idx_interior = 1:length(Model.xVertices_Inner) % loop over all interior
        N_BC_interior = length(xI_InteriorBoundary{idx_interior}(:,1));
        Index_BC_inner{idx_interior+1} = ...
            [(Index_BC_inner{idx_interior}(end)+1):(Index_BC_inner{idx_interior}(end)+N_BC_interior)];
        x1I = [x1I; xI_InteriorBoundary{idx_interior}(:,1)]; 
        x2I = [x2I; xI_InteriorBoundary{idx_interior}(:,2)]; 
    end
    x1I = [x1I; x1I_Interior]; 
    x2I = [x2I; x2I_Interior]; 
    
    Index_BC = 1:N_BC; % can be used by assigning the EBC index
%     Index_BC_inner = N_BC + 1: N_BC + N_BC_Interior;
    
    Index_EBC = Model.CriteriaEBC(x1I_Boundary,x2I_Boundary); % can be used by assigning the EBC index
    Index_NBC = Model.CriteriaNBC(x1I_Boundary,x2I_Boundary);
    Discretization.Index_BC = Index_BC; % for outer
    Discretization.Index_EBC = Index_EBC; % for outer 
    Discretization.Index_NBC = Index_NBC; % for outer
    Discretization.Index_BC_inner = Index_BC_inner;
    Discretization.xVertices_Inner = Model.xVertices_Inner;
else % this is the original code
    Discretization.DomainInteriorExist = 0;
    N_BC = length(x1I_Boundary);
    x1I = [x1I_Boundary; x1I_Interior;]; 
    x2I = [x2I_Boundary; x2I_Interior;];
    Index_BC = 1:N_BC; % can be used by assigning the EBC index
    Index_EBC = Model.CriteriaEBC(x1I_Boundary,x2I_Boundary); % can be used by assigning the EBC index
    Index_NBC = Model.CriteriaNBC(x1I_Boundary,x2I_Boundary);
    Discretization.Index_BC = Index_BC; Discretization.Index_EBC = Index_EBC; Discretization.Index_NBC = Index_NBC;
end

% Creat Nodal List where XP is used for stored the original data
xI = [x1I, x2I];
nP = length(xI);
Discretization.nP = nP;
Discretization.xI = xI;
Discretization.CriteriaEBC = Model.CriteriaEBC;
Discretization.CriteriaNBC = Model.CriteriaNBC;

end
