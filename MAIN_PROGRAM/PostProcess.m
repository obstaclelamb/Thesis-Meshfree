function varargout = PostProcess(RK,Quadrature,Model,Discretization)
% This subroutine visualizes the calculated displacement, strain and stress fields
%  [uhI,Strain,Stress,erroru_atI,errordu_atI]

% load information from discretization
nP = Discretization.nP;
xI = RK.xI; dI = RK.dI; KernelSize = RK.SupportSize;
x1I = xI(:,1); x2I = xI(:,2);
FIGUREPLOSAVE = 1; % true unless one change this to 0

% These comment lines are for save figures, if you want to save figure,
% please un-comment them, as well as the save function following each
% plot below
Model.FolderName = 'RKPM_Solution';
Dir = [Model.FolderName,'\',Quadrature.Integration,'\',Quadrature.Stabilization];
mkdir([Dir]);

% close all exisitng file
close all
%% plot discretization and voronoi cell

switch Quadrature.Integration
  case {'SCNI','DNI'}
    % Plot Voronoi Cell and Local Coordinates
    SHP = Quadrature.SHP; SHPDX1 = Quadrature.SHPDX1; SHPDX2 = Quadrature.SHPDX2;
    VoronoiDiagram = Quadrature.VoronoiDiagram;
    VoronoiCell = VoronoiDiagram.VoronoiCell;
    VerticeCoordinates = VoronoiDiagram.VerticeCoordinates;

    if Model.Plot.Discretization == 1 % if one want to plot discretization
      if FIGUREPLOSAVE % plot the figure and save it
        figure,
        for i = 1:nP %% ith cell hold on
          plot(VerticeCoordinates([VoronoiCell{i}(:);VoronoiCell{i}(1)],1),VerticeCoordinates([VoronoiCell{i}(:);VoronoiCell{i}(1)],2),'k-','linewidth',1.3); hold on
        end
        scatter(xI(:,1),xI(:,2),24,'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]); hold on
        axis([min(xI(:,1)) max(xI(:,1)) min(xI(:,2)) max(xI(:,2))]); pbaspect([abs(min(xI(:,1))-max(xI(:,1))) abs(min(xI(:,2))-max(xI(:,2))) 1])
        xlabel('x1'); ylabel('x2'); title('Nodal Representative Domain');
        set(gca,'FontSize',13);
        saveas(gca,[Dir,'\NodalRepresentativeDomain.jpg']); hold on
      end
    end % end if one Plot.Discretization

  case {'GAUSS'}
    % for Gauss integration, we compute the shape function and nodal
    % points in order to plot it on the nodes
    xQuad = Quadrature.Domain.xQuad;
    xWuqad_BC = Quadrature.BC.xQuad_onBoundary;
    %         SHP = zeros(NP,NP); % Shape Function
    %         SHPDX = zeros(NP,NP); % Deravitive in X
    %         SHPDY = zeros(NP,NP); % Deravitive in Y

    for index_node = 1:nP
      [SHP(index_node,:),SHPDX1(index_node,:),SHPDX2(index_node,:)] = getRKShapeFunction(RK,xI,xI(index_node,:),[1,1,0]);
    end
    if Model.Plot.Discretization == 1 % if one want to plot discretization

      if FIGUREPLOSAVE % plot the figure and save it
        figure, scatter(xI(:,1),xI(:,2),24,'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]); hold on
        %             plot(xI(:,1),xI(:,2),'ko','MarkerSize',5); hold on
        xlabel('x1'); ylabel('x2'); title('Discretization');
        axis([min(xI(:,1)) max(xI(:,1)) min(xI(:,2)) max(xI(:,2))]); pbaspect([abs(min(xI(:,1))-max(xI(:,1))) abs(min(xI(:,2))-max(xI(:,2))) 1]);
        set(gca,'FontSize',13);
        %                 saveas(gca,[Dir,'\Discretization.jpg']); hold on

        figure, plot(xI(:,1),xI(:,2),'ko',xQuad(:,1),xQuad(:,2),'b*',xWuqad_BC(:,1),xWuqad_BC(:,2),'r+','MarkerSize',4); hold on
        xlabel('x1'); ylabel('x2'); title('Background Gauss Points');
        axis([min(xI(:,1)) max(xI(:,1)) min(xI(:,2)) max(xI(:,2))]); pbaspect([abs(min(xI(:,1))-max(xI(:,1))) abs(min(xI(:,2))-max(xI(:,2))) 1])
        legend('Meshfree Points','Background G.I. Points','Boundary G.I. Points','Location','Southeast')
        set(gca,'FontSize',13);
        saveas(gca,[Dir,'\BackgroundGaussCell.jpg']); hold on
      end
    end
end



if FIGUREPLOSAVE % plot the figure and save it
  %% Plot Neighbor Search
  if Model.Plot.Discretization == 1 % if one want to plot discretization

    figure,
    for p_id = fix(1*nP/2)
      % This is just for plotting
      [x_supportedge,y_supportedge] = PlotSupport(xI(p_id,1),xI(p_id,2),KernelSize(p_id),RK.KernelGeometry);
      plot(x_supportedge,y_supportedge,'Color',[0.8,0.8,0.8],'LineWidth',1.5); hold on
    end

    plot(xI(:,1),xI(:,2),'ko','MarkerSize',7); hold on
    % Neighbor search function by rangesearch with different function
    [idx] = rangesearch(xI,xI,KernelSize(fix(1*nP/2)),'Distance',RK.SupportDistanceType);

    plot(xI(idx{fix(1*nP/2)},1),xI(idx{fix(1*nP/2)},2),'ro','MarkerSize',6,'MarkerFaceColor','r'); hold on
    xlabel('x1'); ylabel('x2'); title('Support and Neighbors');
    axis([min(xI(:,1)) max(xI(:,1)) min(xI(:,2)) max(xI(:,2))]); pbaspect([abs(min(xI(:,1))-max(xI(:,1))) abs(min(xI(:,2))-max(xI(:,2))) 1])
    set(get(gca,'XLabel'),'FontSize',16);
    set(get(gca,'YLabel'),'FontSize',16);
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\SupportandNeighbors.jpg']); hold on
  end
end



%% obtain displacement
uh1I = SHP*dI(1:2:end,1);
uh2I = SHP*dI(2:2:end,1);
uhI = [uh1I,uh2I];
% obtain shear stress
% Creat Block Diagonal for C matrix
C = Model.ElasticTensor;
% BLOCKC = kron(eye(nP),C);

% obtain strain and stress field
Strain = zeros(nP,3);
epsilon_11 = SHP*SHPDX1*dI(1:2:end);
epsilon_22 = SHP*SHPDX2*dI(2:2:end);
epsilon_du1dx2 = SHP*SHPDX2*dI(1:2:end);
epsilon_du2dx1 = SHP*SHPDX1*dI(2:2:end);
epsilon_12 = 0.5 * (epsilon_du1dx2 + epsilon_du2dx1);

Strain(:,1) = epsilon_11;
Strain(:,2) = epsilon_22;
Strain(:,3) = 2*epsilon_12;

% Save Stress Field
Stress = zeros(size(Strain));
for idx_nP = 1:nP
  epsilon_11_local = Strain(idx_nP,1);
  epsilon_22_local = Strain(idx_nP,2);
  epsilon_12_local = Strain(idx_nP,3);
  stress_local = C*[epsilon_11_local; epsilon_22_local; epsilon_12_local;];
  Stress(idx_nP,1) = stress_local(1);
  Stress(idx_nP,2) = stress_local(2);
  Stress(idx_nP,3) = stress_local(3);
end
% Stress = BLOCKC*Strain;
sigma_11 = Stress(:,1);
sigma_22 = Stress(:,2);
sigma_12 = Stress(:,3);

if FIGUREPLOSAVE % plot the figure and save it

  % if the domain interior hole exist, for this may version, we always
  % assum there is no hole in the domain, for people who want to try the
  % interior hole in the problem, please check the
  % Pre_GenerateDiscretization.m where we do provide method for
  % generating interior
  if Discretization.DomainInteriorExist
    constraint = []; % constraint for delauny triangulartion
    for idx_interior = 1:length(Model.xVertices_Inner)+1 % loop over all interior
      Index_BC_interior = Discretization.Index_BC_inner{idx_interior}(:);
      constraint = [constraint; ...
        Index_BC_interior,circshift(Index_BC_interior,[-1 0])];
    end
  else
    constraint = [Discretization.Index_BC',circshift(Discretization.Index_BC',[-1 0])];
  end

  %% Using delaunay Triangulation to post processing the scattered data
  tri = delaunayTriangulation(xI(:,1),xI(:,2),constraint);

  if Model.Plot.Displacement == 1
    %% u1
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,uh1I);
    colorbar; colormap(jet); shading interp;
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Displacement, u_{1}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\u1.jpg']); hold on

    %% u2
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,uh2I);
    colorbar; colormap(jet); shading interp; %caxis([up_min up_max])
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Displacement, u_{2}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\u2.jpg']); hold on
  end

  if Model.Plot.Strain == 1
    %% epsilonxx
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,epsilon_11);
    % figure,scatter(XI(:,1)+u1I,XI(:,2)+u2I,36,(u1I),'filled');
    colorbar; colormap(jet); shading interp; % caxis([-0.16 0.16])
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Strain, \epsilon_{11}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\epsilon11.jpg']); hold on

    %% epsilonyy
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,epsilon_22);
    % figure,scatter(XI(:,1)+u1I,XI(:,2)+u2I,36,(u2I),'filled');
    colorbar; colormap(jet); shading interp; %caxis([up_min up_max])
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Strain, \epsilon_{22}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\epsilon22.jpg']); hold on

    %% 2*epsilonxy
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,0.5*epsilon_12);
    colorbar; colormap(jet); shading interp; % caxis([-12E-3 0E-3])
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Strain, \epsilon_{12}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\epsilon12.jpg']); hold on
  end

  if Model.Plot.Stress == 1
    %% sigmaxx
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,sigma_11);
    % figure(3),scatter(XI(:,1)+u1I,XI(:,2)+u2I,36,(u1I),'filled');
    colorbar; colormap(jet); shading interp; % caxis([-0.16 0.16])
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Stress, \sigma_{11}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\sigma11.jpg']); hold on

    %% sigmayy
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,sigma_22);
    % figure(4),scatter(XI(:,1)+u1I,XI(:,2)+u2I,36,(u2I),'filled');
    colorbar; colormap(jet); shading interp; %caxis([up_min up_max])
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Stress, \sigma_{22}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\sigma22.jpg']); hold on

    %% 2*sigmaxy
    figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,0.5*sigma_12);
    colorbar; colormap(jet); shading interp; % caxis([-12E-3 0E-3])
    view(0,90);
    xlabel('x1'); ylabel('x2'); title('Stress, \sigma_{12}');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\sigma12.jpg']); hold on
  end

  if Model.Plot.DeformedConfiguration == 1
    %% deformed configuration and delauny traingle
    figure, s=trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,zeros(size(xI(:,2)))); hold on
    s.EdgeColor = 'k'; s.FaceColor = 'none'; view(0,90);
    scatter(xI(:,1)+uh1I,xI(:,2)+uh2I,24,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]); hold on
    xlabel('x1'); ylabel('x2'); title('Deformed Configuration');
    axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
    pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
    % axis([min(XI(:,1)+u1I) max(XI(:,1)+u1I) min(XI(:,2)+u2I) max(XI(:,2)+u2I)])
    set(gca,'FontSize',13);
    saveas(gca,[Dir,'\DeformedConfiguration.jpg']); hold on
  end

  %% if exact solution exist, compute the nodal error and plot
  if Model.Plot.Error == 1
    if Model.ExactSolution.Exist
      syms x1 x2 n1 n2
      u_exact = Model.ExactSolution.u_exact;
      func_u1exact = matlabFunction(u_exact(1),'Vars',[x1 x2 n1 n2]);
      func_u2exact = matlabFunction(u_exact(2),'Vars',[x1 x2 n1 n2]);
      epsilon_11_ex = diff(u_exact(1),x1);
      epsilon_22_ex = diff(u_exact(2),x2);
      epsilon_12_ex =  0.5*(diff(u_exact(1),x2)+diff(u_exact(2),x1));
      stress = C*[epsilon_11_ex;epsilon_22_ex;2*epsilon_12_ex];
      func_sigma_11 = matlabFunction(stress(1),'Vars',[x1 x2 n1 n2]);
      func_sigma_22 = matlabFunction(stress(2),'Vars',[x1 x2 n1 n2]);
      func_sigma_12 = matlabFunction(stress(3),'Vars',[x1 x2 n1 n2]);
      func_epsilon_11 = matlabFunction(epsilon_11_ex,'Vars',[x1 x2 n1 n2]);
      func_epsilon_22 = matlabFunction(epsilon_22_ex,'Vars',[x1 x2 n1 n2]);
      func_epsilon_12 = matlabFunction(epsilon_12_ex,'Vars',[x1 x2 n1 n2]);
      u1I_exact = func_u1exact(xI(:,1),xI(:,2));
      u2I_exact = func_u2exact(xI(:,1),xI(:,2));
      sigma11_exact = func_sigma_11(xI(:,1),xI(:,2));
      sigma22_exact = func_sigma_22(xI(:,1),xI(:,2));
      sigma12_exact = func_sigma_12(xI(:,1),xI(:,2));
      epsilon11_exact = func_epsilon_11(xI(:,1),xI(:,2));
      epsilon22_exact = func_epsilon_22(xI(:,1),xI(:,2));
      epsilon12_exact = func_epsilon_12(xI(:,1),xI(:,2));
      erroru_atI = sqrt(((uh1I-u1I_exact).^2 + (uh2I-u2I_exact).^2)/...
        sum(((u1I_exact.^2 + (u2I_exact).^2))));
      errordu_atI = sqrt(((epsilon_11-epsilon11_exact).*(sigma_11-sigma11_exact) + ...
        (epsilon_22-epsilon22_exact).*(sigma_22-sigma22_exact) + ...
        (epsilon_12-epsilon12_exact).*(sigma_12-sigma12_exact) )/...
        sum(((epsilon11_exact).*(sigma11_exact) + ...
        (epsilon22_exact).*(sigma22_exact) + ...
        (epsilon12_exact).*(sigma12_exact))));

      figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,erroru_atI);
      colorbar; colormap(jet); shading interp; % caxis([-12E-3 0E-3])
      view(0,90);
      xlabel('x1'); ylabel('x2'); title('Displacement Error');
      axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
      pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
      set(gca,'FontSize',13);
      saveas(gca,[Dir,'\DisplacementError.jpg']); hold on

      figure, trisurf(tri(tri.isInterior(),:), xI(:,1)+uh1I,xI(:,2)+uh2I,errordu_atI);
      colorbar; colormap(jet); shading interp; % caxis([-12E-3 0E-3])
      view(0,90);
      xlabel('x1'); ylabel('x2'); title('Strain Energy Error');
      axis([min(xI(:,1)+uh1I) max(xI(:,1)+uh1I) min(xI(:,2)+uh2I) max(xI(:,2)+uh2I)])
      pbaspect([abs(min(xI(:,1)+uh1I)-max(xI(:,1)+uh1I)) abs(min(xI(:,2)+uh2I)-max(xI(:,2)+uh2I)) 1])
      set(gca,'FontSize',13);
      saveas(gca,[Dir,'\StrainEnergyError.jpg']); hold on

    end
  end
  varargout{1} = uhI;
  varargout{2} = Strain;
  varargout{3} = Stress;
  if nargout == 5
    try
      if exist(erroru_atI, 'var')
        varargout{4} =  erroru_atI;
        varargout{5} = errordu_atI;
      end
    catch
      varargout{4} = 0;
      varargout{5} = 0;
    end
  end
end
end

function [x_edge,y_edge] = PlotSupport(x,y,r,kerneltype)
switch kerneltype
  case 'CIR'
    ang=0:0.01:2*pi;
    xp=r*cos(ang);
    yp=r*sin(ang);
    x_edge = x+xp; y_edge = y+yp;
  case 'REC'
    x_edge = [x-r, x+r, x+r, x-r, x-r];
    y_edge = [y-r, y-r, y+r, y+r, y-r];
end
end