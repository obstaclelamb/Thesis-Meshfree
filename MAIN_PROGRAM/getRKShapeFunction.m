function [PSI, dPSIdx1, dPSIdx2, dPSIdx11, dPSIdx22, dPSIdx12] = getRKShapeFunction(varargin)
% SHAPE FUNCTION OF 2D RK APPROXIMATION
% This subroutine computes the 2D RK shape function and its direct gradients at a given evaluation point
% INPUT PARAMETERS
%    RK:   structure contains all RK information
%    xI:   (RK nodes) vector NP by 2
%    xEva: (Evalutation points): scalar 1 by 2
%    FLAG: [1 1 1]; flag to define which function is looking for. 1: yes
%                   and 0: no
%    Exmaple, if [1 0 0] is used, mean only compute PSI
%    Exmaple, if [1 1 0] is used, mean only compute PSI, dPSIdx1, dPSIdx2
%    Exmaple, if [1 1 1] is used, mean only compute PSI, dPSIdx1, dPSIdx2, dPSIdx1x1, dPSIdx2x2, dPSIdx1x2
%    FLAG_IG: using implicit gradient for the second order gradient


% OUTPUT PARAMETERS
%    PSI     - RK Shape Function
%    dPSIdx1  - First order derivatives of x1
%    dPSIdx2  - First order derivatives of x2
%    Similar for the second order derivatives

% read inputs
RK = varargin{1};
xI = varargin{2};
xEva = varargin{3};
FLAG = [1,1,1]; % [shape function, first derivative, second derivative]
FLAG_IG = 0; % imiplicit gradient is turned off
if length(varargin) == 4
    FLAG = varargin{4}; % flag to use  
elseif length(varargin) == 5
    FLAG = varargin{4};
    FLAG_IG = 1; % flag to use implicity gradient 
end

% because higher derivitive requires lower one so the flag need to be
% adjusted such that they are either [1,1,1] or [1,1,0] or [1,0,0]
if FLAG(3) == 1
    FLAG(1:2) = 1; % flag to use  
elseif FLAG(3) == 0 && FLAG(2) == 1
    FLAG(1) = 1;
end

order  = RK.Order; 
nP = RK.nP; 
type = RK.KernelFunction;
x1I = xI(:,1); x2I = xI(:,2);
x1 = xEva(:,1); x2 = xEva(:,2); 
KernelSize = RK.SupportSize;

% neighbor list Gx construction
switch RK.KernelGeometry
    case 'CIR'
        Gx = transpose(find((bsxfun(@minus,x1I,x1).^2 + bsxfun(@minus,x2I,x2).^2 < KernelSize(:).^2)>0));
    case 'REC'
        Gx = transpose(find((max(bsxfun(@minus,x1I,x1).^2,bsxfun(@minus,x2I,x2).^2) < KernelSize(:).^2)>0));
end
% one can also use the rangesearch function, but after testing, the above
% function is faster than rangesearch in 2D case, 3D may be different
% [NeighborList] = rangesearch([x1I,x2I],[x1,x2],max(KernelSize(:)),'Distance',RK.SupportDistanceType);
% Gx = NeighborList{1}; % neighbor list

d = 2; % two dimensional case

% determine the order of ths problem
switch order
    case {'Constant'}  % constant
        n = 0;
    case {'Linear'} % Linear
        n = 1;
    case {'Quadratic'} % Quadratic
        n = 2;
    otherwise % default is linear
        n = 1;
end

% if FLAG_IG % implicit gradient here is valid for linear basis 
%     n = 1;
% end

% determine numbers of monomial
m = factorial(n+d)/(factorial(n)*factorial(d));

% INITIALIZE H Vectors
H      = zeros (m, nP);  
dHdx1  = zeros (m, nP);
dHdx2  = zeros (m, nP);
dHdx11 = zeros (m, nP);
dHdx22 = zeros (m, nP);
dHdx12 = zeros (m, nP);

% INITIALIZE Kernel functions
PHI      = zeros (1, nP);  
dPHIdx1  = zeros (1, nP);
dPHIdx2  = zeros (1, nP);
dPHIdx11 = zeros (1, nP);
dPHIdx22 = zeros (1, nP);
dPHIdx12 = zeros (1, nP);

% INITIALIZE SHAPE FUNCTION MATRICES
PSI       = zeros(1, nP);
dPSIdx1   = zeros(1, nP);
dPSIdx2   = zeros(1, nP);
dPSIdx11 = zeros(1, nP);
dPSIdx22 = zeros(1, nP);
dPSIdx12 = zeros(1, nP);

% H0 vector
H0 = zeros(m,1); H0(1) = 1;
% calculate H vector, where eps (machine precision) is used to avoid 0^0 = NaN
count = 1;
for k = 0 : n 
    for j = 0 : k 
        i = k - j;
        if FLAG(1) % shape function
        H(count,Gx) = (((x1 - x1I(Gx)+eps)').^i).*(((x2 - x2I(Gx)+eps)').^j);
        end
        if FLAG(2) % first derivative
            if i ~= 0 % in matalb i = 0 will cause 0*inf = NaN when nodal point is used for evaluation shape function
                dHdx1(count,Gx) = (i*((x1 - x1I(Gx)+eps)').^(i-1)).*(((x2 - x2I(Gx)+eps)').^j);
            end
            if j ~= 0
                dHdx2(count,Gx) = (((x1 - x1I(Gx)+eps)').^i).*(j*((x2 - x2I(Gx)+eps)').^(j-1));
            end
        end
        if FLAG(3) % second derivative 
            if i ~= 0 || i ~=1; dHdx11(count,Gx) = (i*(i-1)*((x1 - x1I(Gx)+eps)').^(i-2)).*(((x2 - x2I(Gx)+eps)').^j); end
            if j ~= 0 || j ~=1; dHdx22(count,Gx) = (((x1 - x1I(Gx)+eps)').^i).*(j*(j-1)*((x2 - x2I(Gx)+eps)').^(j-2)); end
            if i ~= 0 || j ~=0; dHdx12(count,Gx) = (i*((x1 - x1I(Gx)+eps)').^(i-1)).*(j*((x2 - x2I(Gx)+eps)').^(j-1)); end
        end
        count = count+1;
    end 
end 


% implicit gradient
if FLAG_IG 
    H1 = zeros(m,1);
    H2 = zeros(m,1);
    if n == 0
    elseif n == 1
        H1(2) = -1;
        H2(3) = -1;
    elseif n == 2
        H1(2) = -1; 
        H2(3) = -1; 
    end
end


% Moment matrix and derivative
Moment   = zeros (m, m);
dMomentdx1= zeros (m, m);
dMomentdx2= zeros (m, m);
dMomentdx11= zeros (m, m);
dMomentdx22= zeros (m, m);
dMomentdx12= zeros (m, m);

% determine the kernel function and the moment matrix. Notice that we can
% do vectorized sum on these functions, but in this form the code is easier
% to manage and do modification
for idx_Gx = 1:length(Gx)
i = Gx(idx_Gx);
[PHI(i), dPHIdx1(i), dPHIdx2(i), dPHIdx11(i), dPHIdx22(i), dPHIdx12(i)] = get_Kernel(type,x1,x2,x1I(i),x2I(i),KernelSize(i),FLAG);
[Moment,dMomentdx1,dMomentdx2,dMomentdx11,dMomentdx22,dMomentdx12] = ...
    getMomentMatrix(i,...
    Moment,dMomentdx1,dMomentdx2,dMomentdx11,dMomentdx22,dMomentdx12,...
    H,dHdx1,dHdx2,dHdx11,dHdx22,dHdx12,...
    PHI,dPHIdx1,dPHIdx2,dPHIdx11,dPHIdx22,dPHIdx12,FLAG);
end


% Calculate the B and correlated matrix
OneVec = ones(m,1);
% the first input is the number of derivative need
% the second to last is A, B, dA, dB, ddA, ddB 
if FLAG(1)
% [B] = getDerivative(0,H,OneVec*PHI);
B = H.*(OneVec*PHI);
end

if FLAG(2)
[dBdx1] = getDerivative(1,H,OneVec*PHI,dHdx1,OneVec*dPHIdx1);
[dBdx2] = getDerivative(1,H,OneVec*PHI,dHdx2,OneVec*dPHIdx2);
end

if (FLAG(3)&&~FLAG_IG)
[dBdx11] = getDerivative(2,H,OneVec*PHI,dHdx1,OneVec*dPHIdx1,dHdx11,OneVec*dPHIdx11);
[dBdx22] = getDerivative(2,H,OneVec*PHI,dHdx2,OneVec*dPHIdx2,dHdx22,OneVec*dPHIdx22);
[dBdx12] = dHdx2 .* (OneVec*dPHIdx1) + H .* (OneVec*dPHIdx12) + dHdx12 .* (OneVec*PHI) + dHdx1.* (OneVec*dPHIdx2);
end

% shape function
    if FLAG(1)
    coef_b  = Moment\H0; % coeficcient b is calculated b = M^{-1}(x)*H0
    PSI(1,:) = coef_b' * B; % shape function                          
    end

% take first derivative
    if FLAG(2)
        if FLAG_IG % imiplicit gradient
            coef_invMH1  = Moment\H1; % coefficient for implicit gradient
            dPSIdx1(1,:) = coef_invMH1' * B; % first order derivatives of shape function to x 

            coef_invMH2  = Moment\H2; % coefficient for implicit gradient
            dPSIdx2(1,:) = coef_invMH2' * B; % first order derivatives of shape function to y  
        else
            % take exact inverse of the moment matrix (M^-1)' = -inv(M)*M'*inv(M);
            dMinvdx1 = -(Moment)\dMomentdx1/(Moment);
            dMinvdx2 = -(Moment)\dMomentdx2/(Moment);

            dcoef_invMdx1H0 = dMinvdx1*H0; % first order derivatives of shape function to x
            dPSIdx1(1,:) = coef_b' * dBdx1 + dcoef_invMdx1H0'*B ;

            dcoef_invMdx2H0 = dMinvdx2*H0; % first order derivatives of shape function to y
            dPSIdx2(1,:) = coef_b' * dBdx2 + dcoef_invMdx2H0'*B ;
        end
    end
    
% take second derivative
    if FLAG(3)
        if FLAG_IG % imiplicit second gradient for take first implicit and second direct direct
            dMinvdx1 = -(Moment)\dMomentdx1/(Moment);
            dMinvdx2 = -(Moment)\dMomentdx2/(Moment);

            dcoef_invMdx1H1 = dMinvdx1*H1; % second order derivatives of shape function to 11
            dPSIdx11(1,:) = coef_invMH1' * dBdx1 + dcoef_invMdx1H1'*B ;

            dcoef_invMdx2H2 = dMinvdx2*H2; % second order derivatives of shape function to 22
            dPSIdx22(1,:) = coef_invMH2' * dBdx2 + dcoef_invMdx2H2'*B ;
            
            dcoef_invMdx2H1 = dMinvdx2*H1; % second order derivatives of shape function to 12
            dPSIdx12(1,:) = coef_invMH1' * dBdx2 + dcoef_invMdx2H1'*B ;
        else
            dMinvdx1x1 = -dMinvdx1*dMomentdx1/(Moment) -(Moment)\dMomentdx11/(Moment) -(Moment)\dMomentdx1*dMinvdx1;
            dMinvdx2x2 = -dMinvdx2*dMomentdx2/(Moment) -(Moment)\dMomentdx22/(Moment) -(Moment)\dMomentdx2*dMinvdx2;
            dMinvdx1x2 = -dMinvdx2*dMomentdx1/(Moment) -(Moment)\dMomentdx12/(Moment) -(Moment)\dMomentdx1*dMinvdx2;

            dcoef_invMdx1H1 = dMinvdx1x1*H0; % second order derivatives of shape function to 11
            dPSIdx11(1,:) = dcoef_invMdx1H0' * dBdx1 + coef_b' * dBdx11 + dcoef_invMdx1H1'*B + dcoef_invMdx1H0'*dBdx1;

            dcoef_invMdx2H2 = dMinvdx2x2*H0; % seoncd order derivatives of shape function to 22
            dPSIdx22(1,:) = dcoef_invMdx2H0' * dBdx2 + coef_b' * dBdx22 + dcoef_invMdx2H2'*B + dcoef_invMdx2H0'*dBdx2 ;

            dcoef_invMdx2H1 = dMinvdx1x2*H0; % seoncd order derivatives of shape function to 12
            dPSIdx12(1,:) = dcoef_invMdx2H0' * dBdx1 + coef_b' * dBdx12 + dcoef_invMdx2H1'*B + dcoef_invMdx1H0'*dBdx2;
            
        end
   
    end
end

%% calculate derivative for genera C = A*B, C,x = A,x*B + A*B,x ... etc
function [C] = getDerivative(varargin)
order_derivative = varargin{1};

    if order_derivative == 0
        A = varargin{2}; B = varargin{3};
        C     = A .* B;
    elseif order_derivative == 1
        A = varargin{2}; B = varargin{3};
        dA = varargin{4}; dB = varargin{5};
        C     = A .* dB + dA .* B;
    elseif order_derivative == 2
        A = varargin{2}; B = varargin{3};
        dA = varargin{4}; dB = varargin{5};
        ddA = varargin{6}; ddB = varargin{7};
        C     = dA .* dB + A .* ddB + ddA .* B + dA .* dB;
    end

end

%% obtain the moment matrix
function [Moment,dMomentdx1,dMomentdx2,dMomentdx11,dMomentdx22,dMomentdx12] = ...
          getMomentMatrix(j,...
          Moment,dMomentdx1,dMomentdx2,dMomentdx11,dMomentdx22,dMomentdx12,...
          H,dHdx1,dHdx2,dHdx11,dHdx22,dHdx12,...
          PHI,dPHIdx1,dPHIdx2,dPHIdx11,dPHIdx2x2,dPHIdx12,...
          FLAG)
% moment
if FLAG(1)
Moment    = Moment    + PHI(j) * H(:,j) * H(:,j)';
end

% first derivative
if FLAG(2)
dMomentdx1 = dMomentdx1 + dPHIdx1(j) * H(:,j) * H(:,j)' + PHI(j) * dHdx1(:,j) *H(:,j)' + PHI(j) * H(:,j)*dHdx1(:,j)';
dMomentdx2 = dMomentdx2 + dPHIdx2(j) * H(:,j) * H(:,j)' + PHI(j) * dHdx2(:,j) *H(:,j)' + PHI(j) * H(:,j)*dHdx2(:,j)';
end

% second derivative
if FLAG(3)
dMomentdx11 = dMomentdx11 + dPHIdx11(j) * H(:,j) * H(:,j)' + dPHIdx1(j) * dHdx1(:,j) * H(:,j)' + dPHIdx1(j) * H(:,j) * dHdx1(:,j)' + ...
                        dPHIdx1(j) * dHdx1(:,j) *H(:,j)' + PHI(j) * dHdx11(:,j) *H(:,j)' + PHI(j) * dHdx1(:,j) *dHdx1(:,j)' + ...
                        dPHIdx1(j) * H(:,j)*dHdx1(:,j)' + PHI(j) * dHdx1(:,j)*dHdx1(:,j)' + PHI(j) * H(:,j)*dHdx11(:,j)';
dMomentdx22 = dMomentdx22 + dPHIdx2x2(j) * H(:,j) * H(:,j)' + dPHIdx2(j) * dHdx2(:,j) * H(:,j)' + dPHIdx2(j) * H(:,j) * dHdx2(:,j)' + ...
                        dPHIdx2(j) * dHdx2(:,j) *H(:,j)' + PHI(j) * dHdx22(:,j) *H(:,j)' + PHI(j) * dHdx2(:,j) *dHdx2(:,j)' +...
                        dPHIdx2(j) * H(:,j)*dHdx2(:,j)' + PHI(j) * dHdx2(:,j)*dHdx2(:,j)' + PHI(j) * H(:,j)*dHdx22(:,j)';
dMomentdx12 = dMomentdx12 + dPHIdx12(j) * H(:,j) * H(:,j)' + dPHIdx1(j) * dHdx2(:,j) * H(:,j)' + dPHIdx1(j) * H(:,j) * dHdx2(:,j)'...
                      + dPHIdx2(j) * dHdx1(:,j) *H(:,j)' + PHI(j) * dHdx12(:,j) *H(:,j)' + PHI(j) * dHdx1(:,j) *dHdx2(:,j)'...
                      + dPHIdx2(j) * H(:,j)*dHdx1(:,j)' + PHI(j) * dHdx2(:,j)*dHdx1(:,j)' + PHI(j) * H(:,j)*dHdx12(:,j)';
end

end

%% obtain kernel function
function [w, dwdx1, dwdx2, dwdx11, dwdx22, dwdx12] = get_Kernel(type,x1,x2,x1I,x2I,a,FLAG)
% for saving space, w is expression for phi
w=0;
dwdz=0;
dwdzz=0;
dzdx1=0;
dzdx2=0;
dzdx11=0;
dzdx22=0;
dzdx12=0;

% evalutati kernel distance between x and xI
z  = sqrt((x1-x1I).^2+(x2-x2I).^2)/a;
% nominator for z, eps is to prevent singularity
z_nominator = (a.^2*z)+eps;
% first derivative
if FLAG(2) 
dzdx1 = (x1-x1I)/(z_nominator);
dzdx2 = (x2-x2I)/(z_nominator);
end
% second derivative
if FLAG(3)
dzdx11 = (1*z_nominator - (x1-x1I).*a.^2*dzdx1)/(z_nominator)^2;
dzdx22 = (1*z_nominator - (x2-x2I).*a.^2*dzdx2)/(z_nominator)^2;
dzdx12 = (-(x1-x1I).*a.^2*dzdx2)/(z_nominator)^2;
end 
% evaluate kernel function with it's derivative respect to z
switch type
    case 'HVSIDE'
        [w,dwdz,dwdzz] = HVSIDE(z);
    case 'SPLIN1'
        [w,dwdz,dwdzz] = SPLIN1(z);
    case 'SPLIN2'
        [w,dwdz,dwdzz] = SPLIN2(z);
    case 'SPLIN3'
        [w,dwdz,dwdzz] = SPLIN3(z);
    case 'SPLIN4'
        [w,dwdz,dwdzz] = SPLIN4(z);
    case 'SPLIN5'
        [w,dwdz,dwdzz] = SPLIN5(z);
    otherwise
        [w,dwdz,dwdzz] = SPLIN3(z);
end
% obtain kernel function
dwdx1 = dwdz .* dzdx1;
dwdx2 = dwdz .* dzdx2;
dwdx11 = dwdz .* dzdx11 + dwdzz .* dzdx1.* dzdx1;
dwdx22 = dwdz .* dzdx22 + dwdzz .* dzdx2.* dzdx2 ;
dwdx12 = dwdz .* dzdx12 + dwdzz .* dzdx1.* dzdx2;

end

%% Some Kernel function
function [w,dwdr,dwdrr] = HVSIDE(z)
if (z>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr  = 0.0;
else
	w     = 1;
	dwdr  = 0;
    dwdrr  = 0;
end
end

function [w,dwdr,dwdrr] = SPLIN1(z)
if (z>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr  = 0.0;
else
	w     = 1-z;
	dwdr  = -1;
    dwdrr  = 0;
end
end

function [w,dwdr,dwdrr] = SPLIN2(z)
if (z>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr = 0;
elseif (z>1/3)
	w     = 3/2 - 3*z + 3/2*z^2;
	dwdr  = -3 + 3*z;
    dwdrr = 3*z;
else
	w     = 1 - 3*z^2;
	dwdr  = -6*z; 
    dwdrr = -6;
end
end

function [w,dwdr,dwdrr] = SPLIN3(z)
if (z>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr  = 0.0;
elseif (z<=0.5)
   w     = 2/3 - 4*z^2 + 4*z^3;
   dwdr  = -8*z + 12*z^2;
   dwdrr  = -8 + 24*z;
else
   w     = 4/3 - 4*z + 4*z^2 - 4*z^3/3;
   dwdr  = -4 + 8*z -4*z^2;
   dwdrr  = 8 - 8*z;
end
end

function [w,dwdr,dwdrr] = SPLIN4(z)
if (z>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr = 0;
elseif (z > 3.0/5.0)
	w     = 125.0/46.0 - 250.0/23.0*z + 375.0/23.0*z^2 - ...
            250.0/23.0*z^3 + 125.0/46.0*z^4;
	dwdr  = - 250.0/23.0 + 2*375.0/23.0*z - ...
            3*250.0/23.0*z^2 + 4*125.0/46.0*z^3;
    dwdrr = 2*375.0/23.0 - 2*3*250.0/23.0*z + 3*4*125.0/46.0*z^2;
elseif (z > 1.0/5.0) && (z <= 3.0/5.0)
	w     = 22.0/23.0 + 20.0/23.0*z - 300.0/23.0*z^2 + ...
            500.0/23.0*z^3 - 250.0/23.0*z^4;
	dwdr  = + 20.0/23.0 - 2*300.0/23.0*z + ...
            3*500.0/23.0*z^2 - 4*250.0/23.0*z^3; 
    dwdrr  = 2*300.0/23.0 + 2*3*500.0/23.0*z - 3*4*250.0/23.0*z^2; 
else
	w     = 1.0 - 150.0/23.0*z^2 + 375.0/23.0*z^4;
	dwdr  = - 2*150.0/23.0*z + 4*375.0/23.0*z^3;
    dwdrr  = - 2*150.0/23.0 + 3*4*375.0/23.0*z^2;
end
end

function [w,dwdr,dwdrr] = SPLIN5(z)
if (z>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr = 0.0;
elseif (z > 2.0/3.0)
	w     = 81.0/22.0 - 405.0/22.0*z + 405.0/11.0*z^2 - ...
            405.0/11.0*z^3 + 405.0/22.0*z^4 - ...
            81.0/22.0*z^5;
	dwdr  = - 405.0/22.0 + 2*405.0/11.0*z - ...
            3*405.0/11.0*z^2 + 4*405.0/22.0*z^3 - ...
            5*81.0/22.0*z^4;
    dwdrr = 2*405.0/11.0 - 2*3*405.0/11.0*z + ...
            3*4*405.0/22.0*z^2 - 4*5*81.0/22.0*z^3;
elseif (z > 1.0/3.0) && (z <= 2.0/3.0)
	w     = 17.0/22.0 + 75.0/22.0*z - 315.0/11.0*z^2 + ...
            675.0/11.0*z^3 - 1215.0/22.0*z^4 + ...
            405.0/22.0*z^5;
	dwdr  = + 75.0/22.0 - 2*315.0/11.0*z + ...
            3*675.0/11.0*z^2 - 4*1215.0/22.0*z^3 + ...
            5*405.0/22.0*z^4;
    dwdrr  = - 2*315.0/11.0 + ...
            2*3*675.0/11.0*z - 3*4*1215.0/22.0*z^2 + ...
            4*5*405.0/22.0*z^3;    
else
	w     = 1.0 - 90.0/11.0*z^2 + 405.0/11.0*z^4 - ...
            405.0/11.0*z^5;
	dwdr  = - 2*90.0/11.0*z + 4*405.0/11.0*z^3 - ...
            5*405.0/11.0*z^4;
    dwdrr  = - 2*90.0/11.0 + 3*4*405.0/11.0*z^2 - ...
            4*5*405.0/11.0*z^3;
end
end