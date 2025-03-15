function [HXS,DHXS,DDHXS,H0,DH0,DDH0]=get_H_matrix(XS,DEG)
% get vector of basis function
np = length(XS);
x0 = zeros(1,np);
x1 = ones(1,np);
H0 = [x1;zeros(DEG,np)];      %get H(0)
DH0 = zeros(DEG+1,np);        %get dH(0)
DDH0 = zeros(DEG+1,np);       %get ddH(0)

n0 = (0:DEG);% 0123
n1 = (1:DEG);%  123
HXS   = XS.^(n0');            %get H(X-XI)
DHXS  = [x0;HXS(n1,:).*(n1')];  %get dH(X-XI)
DDHXS = [x0;DHXS(n1,:).*(n1')]; %get ddH(X-XI)
%{ 
%% slow
if DEG<=0
  error('degree error in GETHMAT');
end
DEG = floor(DEG);   %% in case of float value
syms x;
HXS = x.^(0:DEG)';
DHXS = diff(HXS,x);
DDHXS = diff(DHXS,x);
HXS = subs(HXS,x,XS);
DHXS = subs(DHXS,x,XS);
DDHXS = subs(DDHXS,x,XS);

np = length(XS);
H0 = [ones(1,np);zeros(DEG,np)];
DH0 = zeros(DEG+1,np);
DDH0 = zeros(DEG+1,np);
%%  deprecated
HXS  =[x1;    XS;    XS.^2;    XS.^3;     XS.^4;     XS.^5;     XS.^6];
DHXS =[x0;    x1; 2.*XS   ; 3.*XS.^2;  4.*XS.^3;  5.*XS.^4;  6.*XS.^5];
DDHXS=[x0;    x0; 2.*X1   ; 6.*XS   ; 12.*XS.^2; 20.*XS.^3; 30.*XS.^4];

HXS = HXS(1:DEG+1,:);
DHXS = DHXS(1:DEG+1,:);
DDHXS = DDHXS(1:DEG+1,:);
%% original
if DEG==0;
    HXS=[1.]';
    DHXS=[0]';
    DDHXS=[0]';
    H0=[1.]';
    DH0=[0]';
    DDH0=[0]';
elseif DEG==1;
    HXS=[1.,XS]';
    DHXS=[0,1.]';
    DDHXS=[0,0]';
    H0=[1.,0]';
    DH0=[0,0]';
    DDH0=[0,0]';
elseif DEG==2;
    HXS=[1.,XS,XS^2]';
    DHXS=[0,1.,2.*XS]';
    DDHXS=[0,0,2.]';
    H0=[1.,0,0]';
    DH0=[0,0,0]';
    DDH0=[0,0,0]';
elseif DEG==3;
    HXS=[1.,XS,XS^2,XS^3]';
    DHXS=[0,1.,2*XS,3.*XS^2]';
    DDHXS=[0,0,2.,6.*XS]';
    H0=[1.,0,0,0]';
    DH0=[0,0,0,0]';
    DDH0=[0,0,0,0]';
elseif DEG==4;
    HXS=[1.,XS,XS^2,XS^3,XS^4]';
    DHXS=[0,1.,2.*XS,3.*XS^2,4.*XS^3]';
    DDHXS=[0,0,2.,6.*XS,12.*XS^2]';
    H0=[1,0,0,0,0]';
    DH0=[0,0,0,0,0]';
    DDH0=[0,0,0,0,0]';
elseif DEG==5;
    HXS=[1.,XS,XS^2,XS^3,XS^4,XS^5]';
    DHXS=[0,1.,2.*XS,3.*XS^2,4.*XS^3,5.*XS^4]';
    DDHXS=[0,0,2.,6.*XS,12.*XS^2,20.*XS^3]';
    H0=[1,0,0,0,0,0]';
    DH0=[0,0,0,0,0,0]';
    DDH0=[0,0,0,0,0,0]';
elseif DEG==6;
    HXS=[1.,XS,XS^2,XS^3,XS^4,XS^5,XS^6]';
    DHXS=[0,1.,2.*XS,3.*XS^2,4.*XS^3,5.*XS^4,6.*XS^5]';
    DDHXS=[0,0,2.,6.*XS,12.*XS^2,20.*XS^3,30.*XS^4]';
    H0=[1,0,0,0,0,0,0]';
    DH0=[0,0,0,0,0,0,0]';
    DDH0=[0,0,0,0,0,0,0]';
else
    error ('please choose DEG=(0-6)\n')
end

%}