function [w,x]=getGaussPt(N,varargin)
% Legendre-Gauss Quadrature. 
% Generate Gauss Quadrature point for interval [min, max]
% the default interval is set to [-1,1]
% GetGaussPt(order)
% GetGaussPt(order, [min, max])
% GetGaussPt(order, [min], [max])

% modified from: 
% https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
% Written by Aiyung Chen - 2024/11/25
switch length(varargin) 
  case 0
    a = -1; b = 1;
  case 1
    if any(size(varargin{1}) ~= [1,2])
      a = -1; b = 1;
      warning('Second parameter of GetGaussPt not in [min,max] format')
    else
      a = varargin{1}(1); b = varargin{1}(2);
    end
  case 2
    a = varargin{1}; b = varargin{2};
  otherwise
    error('Parameter error in lgwt');
end

N=N-1; N1=N+1; N2=N+2;
xu=linspace(-1,1,N1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);
% Derivative of LGVM
Lp=zeros(N1,N2);
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0 = 2;
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    L(:,1)=1;
    Lp(:,1)=0;
    L(:,2)=y;
    Lp(:,2)=1;
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end
% Linear map from[-1,1] to [a,b]
x=flip((a*(1-y)+b*(1+y))/2)';  
% Compute the weights
w=flip((b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2)';