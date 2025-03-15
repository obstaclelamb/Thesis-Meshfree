% MATLAB code for meshfree method project 1
% Comment/uncomment specific question line for demostration
% for more info type:
% 
% >>help RK1D
% 
% Aiyung Chen, 2024/11/30
close all;clear;clc;
f = @(x)(ones(size(x))); f_str = 'constant'; Q = 'a';
% f = @(x)(x.*x);          f_str = 'x^{2}'; 
% f = @(x)(cos(x));        f_str = 'cos(x)';Q = 'b';
% f = @(x)(x.^3);          f_str = 'x^{3}';Q = 'c';
% f = @(x)(exp(2.*x));     f_str = 'e^{2x}';Q = 'd';

% node = sort(rand(1,11)); node =(node-node(1))./(node(end)-node(1));
node = linspace(0,1,11); st = 'Structured Node';
% node = [0,0.1,0.18,0.24,0.33,0.45,0.61,0.65,0.8,0.86,1.0];st = 'Non-structured Node';
% node = linspace(0,2*pi,11);
dia = 3;  % dilation
DEG = 1;  % degree of basis vector H 
% kernel = 'Tent';
% kernel = 'Heaviside';
% kernel = 'Quadratic B';
% kernel = 'Cubic B';
% kernel = 'Quartic B'; X
% kernel = 'Quartic';
kernel = 'Quintic B';
% kernel = 'Gaussian';

f_RK = RK1D(f,node,dia,DEG,kernel);

%% partition of unity
% f_RK.shape_fun = zeros(201,11);
% tent = [linspace(0,1,21)';linspace(1,0,21)'];
% tent(21)=[];
% 
% for i=2:10
%   f_RK.shape_fun(20*(i-2)+1:20*(i)+1,i)=tent;
% end
% f_RK.shape_fun(1:21,1)=linspace(1,0,21)';
% % f_RK.shape_fun(1:21,1)=linspace(0,1,21)';
% f_RK.shape_fun(181:201,11)=linspace(0,1,21)';
% f_RK.X = linspace(0,1,201);
%% plot the result
X = f_RK.X;
str = [', DEG: ',num2str(DEG),', dia: ',num2str(dia),', kernel: ',kernel];
% plot shape function and nodes
plotShapeFunction(['$\Psi$',str],f_RK.shape_fun,X);xlabel(st);
scatter(node,zeros(size(node)),50,'r','filled');
scatter(X,sum(f_RK.shape_fun,2)');
ylim([-0.2,1.3]);
% legend(lengend_text(node));  % make legend of each node invisible
% plot derivatives of shape function
plotShapeFunction(['$\frac{d}{dx}\Psi$',str],f_RK.Dshape_fun,X);xlabel(st);
plotShapeFunction(['$\frac{d^2}{dx^2}\Psi$',str],f_RK.DDshape_fun,X);xlabel(st);
% plot completeness
plotCompleteness(f_str,f_RK)

%% output image
% % % pngnamess = replace(char(datetime("now")),":","-");
% name = [Q,',DEG-',num2str(DEG),',',kernel,',s'];
% for i = 1:3
%   figure(i);
% %   pngmaker(string([name,'--',num2str(i),' ',pngnamess]));
%   pngmaker(string([name,'--',num2str(i)]));
% end
% close all;
% %% other
% function T = lengend_text(node)
% np = length(node);
% T = cell(1,np+2);
% T(end-1)={'node'};
% T(end) = {'unity'};
%   for i=1:np
%     T(i) = {''};
%   end
% end