function [N0,N1,N2] = getMLS_ShapeFun(x,DEG,xv,a,spline_idx)
% x:      single val,  point to be interpolated
% DEG     single val,  degree of shape function
% xv      array,       node locatoin
% a       array,       dilation radius of each node (xv)
% spline  single val,  the name of spline we are using
% note: this routine is not optimized but could be
%%
np = length(xv);
%initalize shape function vectors
N0=zeros(1,np); N1=zeros(1,np); N2=zeros(1,np);
%initalize moment matrix and its derivatives
M=zeros((DEG+1),(DEG+1)); dM=M; ddM=M;
IM=eye(DEG+1);

% difference of one node(x) to other node(xv)
XS = x - xv;
ZL=abs(XS)./a;
resolu = 1.0e-15;
DZDX = ((ZL>resolu).*sign(XS) + (ZL<=resolu).*1)./a;
[mHXS,mDHXS,mDDHXS,mH0,mDH0,mDDH0]=get_H_matrix(XS,DEG); 
[PHIZL,DPHIZL,DDPHIZL]=getPhiDD(ZL,spline_idx);  
phi_arr=PHIZL./a;
dphi_arr=DPHIZL.*DZDX./a;
ddphi_arr=(DDPHIZL.*DZDX.^2)./a;
%first loop: get moment matrix
for j = 1:np
  phi = phi_arr(j);    dphi = dphi_arr(j);    ddphi = ddphi_arr(j);
  HXS = mHXS(:,j);    DHXS = mDHXS(:,j);    DDHXS = mDDHXS(:,j);

  M  =  M+   HXS*HXS'*phi;

  dM = dM+  DHXS*HXS'*phi+HXS*DHXS'*phi+HXS*HXS'*dphi;

  ddM=ddM+ DDHXS*HXS'*phi+DHXS*DHXS'*phi+DHXS*HXS'*dphi+ ...
    DHXS*DHXS'*phi+HXS*DDHXS'*phi+HXS*DHXS'*dphi+ ...
    DHXS*HXS'*dphi+HXS*DHXS'*dphi+HXS*HXS'*ddphi;
end
%get derivatives of the inverses of moment matricies
invm=IM/M;
dinvm=-invm*dM*invm;
ddinvm=-invm*ddM*invm-invm*dM*dinvm-invm*dM*dinvm;
%second loop: get shape functions
for j=1:np
  phi = phi_arr(j);    dphi = dphi_arr(j);    ddphi = ddphi_arr(j);
  HXS = mHXS(:,j);    DHXS = mDHXS(:,j);    DDHXS = mDDHXS(:,j);
  H0 = mH0(:,j);    DH0 = mDH0(:,j);    DDH0 = mDDH0(:,j);
  N0(j)=H0'*invm*HXS*phi;

  N1(j)= DH0'* invm* HXS* phi+ H0'* dinvm* HXS* phi+ H0'* invm* DHXS* phi+ H0'* invm* HXS* dphi;

  N2(j)=DDH0'* invm* HXS* phi+DH0'* dinvm* HXS* phi+DH0'* invm* DHXS* phi+DH0'* invm* HXS* dphi+ ...
         DH0'*dinvm* HXS* phi+ H0'*ddinvm* HXS* phi+ H0'*dinvm* DHXS* phi+ H0'*dinvm* HXS* dphi+ ...
         DH0'* invm*DHXS* phi+ H0'* dinvm*DHXS* phi+ H0'* invm*DDHXS* phi+ H0'* invm*DHXS* dphi+ ...
         DH0'* invm* HXS*dphi+ H0'* dinvm* HXS*dphi+ H0'* invm* DHXS*dphi+ H0'* invm* HXS*ddphi;
end