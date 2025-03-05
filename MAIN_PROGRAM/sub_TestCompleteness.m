function [PASS_ORNOT] = sub_TestCompleteness(ORDER,XK,XEVA,SHP,SHPDX,SHPDY,drawplot)
% This subroutine tests the reproducing condition of RK shape functions and gradients.

% Here is an Example of testing the partition of unity and nullity
% for approximate Sin Cos function
UI = sin(XK(:,1)).*cos(XK(:,2));
dUIdX = cos(XK(:,1)).*cos(XK(:,2)); 
dUIdY = -sin(XK(:,1)).*sin(XK(:,2));

    switch drawplot
        case 'ON'
            % partition of uSHPty
            figure(6),scatter3(XK(:,1),XK(:,2),UI,36,'k','filled'); hold on
            figure(6),scatter3(XK(:,1),XK(:,2),SHP*UI,36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('U'); title('U vs U^{h}');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
%             saveas(gca,'U.jpg'); hold on

            figure(7),scatter3(XK(:,1),XK(:,2),dUIdX,36,'k','filled'); hold on
            figure(7),scatter3(XK(:,1),XK(:,2),SHPDX*UI,36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('dUdx'); title('U_{,x} vs U_{,x}^{h}');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
%             saveas(gca,'dUdX.jpg'); hold on

            figure(8),scatter3(XK(:,1),XK(:,2),dUIdY,36,'k','filled'); hold on
            figure(8),scatter3(XK(:,1),XK(:,2),SHPDY*UI,36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('dUdy'); title('U_{,y} vs U_{,y}^{h}');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
%             saveas(gca,'dUdY.jpg'); hold on
            
            %% Check reproducing condition
            figure(9),scatter3(XK(:,1),XK(:,2),sum(SHP')-1,36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('\Sigma\Phi_{I}-1'); title('Constant Consistency');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
%             saveas(gca,'ConstantConsistency.jpg'); hold on

            figure(10),scatter3(XK(:,1),XK(:,2),SHP*XK(:,1)-XEVA(:,1),36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('\Sigma\Phi_{I}X_{I}-X'); title('Linear Consistency');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
%             saveas(gca,'LinearConsistency.jpg'); hold on

            figure(11),scatter3(XK(:,1),XK(:,2),SHPDX*XK(:,1)-1,36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('\Sigma\Phi_{I,x}X_{I}-1'); title('Differential Consistency');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
%             saveas(gca,'DifferentialConsistency.jpg'); hold on
    
            figure(12),scatter3(XK(:,1),XK(:,2),SHPDY*XK(:,2)-1,36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('\Sigma\Phi_{I,y}X_{I}-1'); title('Differential Consistency');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
            
            figure(13),scatter3(XK(:,1),XK(:,2),SHPDY*XK(:,1),36,'ro'); hold on
            xlabel('x'); ylabel('y'); zlabel('\Sigma\Phi_{I,x}X_{I}-1'); title('Differential Consistency');
            set(get(gca,'XLabel'),'FontSize',16);
            set(get(gca,'YLabel'),'FontSize',16);
            set(gca,'FontSize',13);
            
%             saveas(gca,'DifferentialConsistency.jpg'); hold on
            
        otherwise
    end
    
    switch ORDER
        case 'Constant'
            if (norm(sum(SHP')-1) < 1E-7) 
                PASS_ORNOT = 1;
                if (sum(sum(isnan(SHPDX)))>0) || (sum(sum(isnan(SHPDY)))>0) || (sum(sum(isnan(SHP)))>0)
                    PASS_ORNOT = 0;
                end
            else
                PASS_ORNOT = 0;
            end
        case {'Linear','Quadratic'}
            if (norm(SHPDX*XK(:,1)-1) < 1E-7) && (norm(SHP*XK(:,1)-XEVA(:,1)) < 1E-7) && (norm(sum(SHP')-1) < 1E-7) 
                PASS_ORNOT = 1;
                if (sum(sum(isnan(SHPDX)))>0) || (sum(sum(isnan(SHPDY)))>0) || (sum(sum(isnan(SHP)))>0)
                    PASS_ORNOT = 0;
                end
            else
                PASS_ORNOT = 0;
            end
    end
    if ~PASS_ORNOT
    disp('The discretization will make the RK shape function violate the completeness.');
    disp('Please check numerical setting, eg. support size.');
    disp(['Constant Completeness: ',num2str(norm(sum(SHP')-1))]);
    disp(['Linear Completeness: ',num2str(norm(SHP*XK(:,1)-XEVA(:,1)))]);
    disp(['Differential Completeness: ',num2str(norm(SHPDX*XK(:,1)-1))]);
    end
end

