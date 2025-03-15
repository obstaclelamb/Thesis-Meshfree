function result = ...
  RK1D(exact_function, nodes, dilation, basis_order, spline_idx)
% Demostration of Project 1.
% >> RK1D(exact_function, node, dilation, basis_order, spline_idx)
%  exact fun:    function handle of the one to be reproduced.
%  node:         nodes array. can be structured or non-structured
%  dilation:     how many nodes contained by support (c in a = ch)
%  basis_order:  literally
%  spline_idx:   one of the following:
%        Tag           Continuity
%       'Tent'           C0
%       'Heaviside'      C -1
%       'Quadratic B'    C2
%       'Cubic B'        C3
%       'Quartic B'      C1
%       'Quartic'        C2
%       'Quintic B'      C4
%       'Gaussian'       C infty
%
%  domain is considered [node(1), node(end)] by default
% 
% Example input:
%  r = RK1D(@(x)(exp(x)), linspace(0,1,11), 2.5, 4, 'Quintic B'); 

% Aiyung chen, 2024/11/30
f = exact_function;
DEG = basis_order;
dx_min = min(diff(nodes));
dx_max = max(diff(nodes));
support_size = ones(size(nodes))*dilation*dx_max;
samplePt_x = nodes(1):dx_min/20:nodes(end);
C = zeros(1,length(samplePt_x));
DC = C; DDC = C;
N0_BIG = zeros(length(samplePt_x),length(nodes));
N1_BIG = zeros(length(samplePt_x),length(nodes));
N2_BIG = zeros(length(samplePt_x),length(nodes));
for i=1:length(samplePt_x)
    xp = samplePt_x(i);
    [N0,N1,N2] = getMLS_ShapeFun(xp,DEG,nodes,support_size,spline_idx);
    C(i) = sum(N0.*f(nodes));
    DC(i) = sum(N1.*f(nodes));
    DDC(i) = sum(N2.*f(nodes));
    %plot the shape functions
    N0_BIG(i,:) = N0;
    N1_BIG(i,:) = N1;
    N2_BIG(i,:) = N2; 
end
%% return
result.X = samplePt_x;
result.f = exact_function;
result.DEG = DEG;
result.kernel = spline_idx;
result.shape_fun   = N0_BIG;
result.Dshape_fun  = N1_BIG;
result.DDshape_fun = N2_BIG;
result.C = C;
result.DC = DC;
result.DDC = DDC;
end