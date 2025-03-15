function [PHIZL,DPHIZL,DDPHIZL] = getPhiDD(ZL,spline_tag)
% return kernel function, input normalized deviation and spline tag
switch(spline_tag)
  case 'Tent'        % C0 continuety
    PHIZL=  (ZL<=1).*(1.0-ZL);
    DPHIZL= (ZL<=1).*(-1.0);
    DDPHIZL= zeros(size(ZL));
  case 'Heaviside'   % C-1 continuety
    PHIZL=  (ZL<=1).*(1.0);
    DPHIZL=  zeros(size(ZL));
    DDPHIZL= zeros(size(ZL));
  case 'Quadratic B' % C2
    ZL=ZL.*1.5;
    condi1 = ZL<=0.5;
    condi2 = not(condi1)&(ZL<=1.5);
    PHIZL =   condi1.*(3.0/4.0-ZL.^2) + condi2.*(9.0/8.0-3.0/2.0.*ZL+0.5.*ZL.^2);
    DPHIZL =  condi1.*(-2.0.*ZL)      + condi2.*(-3.0/2.0+ZL);
    DDPHIZL = condi1.*(-2.0)          + condi2.*(1.0);
  case 'Cubic B'    % C3 continuety
    condi1 = ZL<=0.5;
    condi2 = not(condi1)&(ZL<=1.0);
    PHIZL=   condi1.*( 2.0/3.0 -4.0.*ZL.^2 +4.0.*ZL.^3) ...
      + condi2.*(4.0/3.0 -4.0.*ZL +  4.0.*ZL.^2  - 4.0/3.0.*ZL.^3);
    DPHIZL=  condi1.*(-8.0.*ZL+12.0.*ZL.^2) + condi2.*(-4.0 +8.0.*ZL -4.0.*ZL.^2);
    DDPHIZL= condi1.*(-8.0 + 24.0.*ZL)      + condi2.*(8.0 -8.0.*ZL);
  case 'Quartic B'
    error('doesn''t have time');
  case 'Quartic'    % C2 continuety
    condi = (ZL<=1);
    PHIZL =   condi.*(1.0 -6.0.*ZL.^2 +8.0.*ZL.^3 -3.0.*ZL.^4);
    DPHIZL =  condi.*(-12.0.*ZL+ 24.0.*ZL.^2 -12.0*ZL.^3);
    DDPHIZL = condi.*(-12.0   + 48.0.*ZL -36.0.*ZL.^2);
  case 'Quintic B'    % (11/20) * quintic B-spline (which is not C4)
    condi1=(ZL<=1.0/3.0);
    condi2=not(condi1)&(ZL<=2.0/3.0);
    condi3=not(condi1|condi2)&(ZL<=1.0);
    PHIZL=condi1.*(11.0/20.0-9.0/2.0.*ZL.^2+81.0/4.0.*ZL.^4-81.0/4.0.*ZL.^5) ...
      + condi2.*(17.0/40.0+15.0/8.0.*ZL-63.0/4.0.*ZL.^2+135.0/4.0.*ZL.^3-243.0/8.0.*ZL.^4+81.0/8.0.*ZL.^5) ...
      + condi3.*(81.0/40.0-81.0/8.0.*ZL+81.0/4.0.*ZL.^2-81.0/4.0.*ZL.^3+81.0/8.0.*ZL.^4-81.0/40.0.*ZL.^5);
    DPHIZL=condi1.*(-9.0.*ZL+81.0.*ZL.^3-5.0.*81.0/4.0.*ZL.^4) ...
      + condi2.*(15.0/8.0-63.0/2.0.*ZL+3.0.*135.0/4.0.*ZL.^2-243.0/2.0.*ZL.^3+5.0.*81.0/8.0.*ZL.^4) ...
      +condi3.*(-81.0/8.0+81.0/2.0.*ZL-3.0.*81.0/4.0.*ZL.^2+81.0/2.0.*ZL.^3-5.0.*81.0/40.0.*ZL.^4);
    DDPHIZL=condi1.*(-9.0+3.0.*81.0.*ZL.^2-5.0.*81.0.*ZL.^3)...
      +condi2.*(-63.0/2.0+3.0.*135.0/2.0.*ZL-3.0.*243.0/2.0.*ZL.^2+5.0.*81.0/2.0.*ZL.^3)...
      +condi3.*(+81.0/2.0-3.0.*81.0/2.0.*ZL+3.0.*81.0/2.0.*ZL.^2-5.0.*81.0/10.0.*ZL.^3);
  case 'Gaussian'    % C0 continuety at edges, C infinity in support
    cg=6.0;
    condi = (ZL<=1);
    PHIZL=condi.*exp(-cg.*ZL);
    DPHIZL=condi.*(-cg.*PHIZL);
    DDPHIZL=condi.*(-cg.*DPHIZL);
    PHIZL0 = exp(-cg*1);
    PHIZL=PHIZL-PHIZL0;
  otherwise
    error('u:stuffed:it', ...
      ['error in getPhiDD function:\n', ...
      '  getPhiDD(spline_tag)\n', ...
      'The kernel function index should be one of following:\n', ...
      ' Tag           Continuity\n', ...
      '''Tent''           C0\n', ...
      '''Heaviside''      C1\n', ...
      '''Quadratic B''    C2\n', ...
      '''Cubic B''        C3\n', ...
      '''Quartic B''      C1\n', ...
      '''Quartic''        C2\n', ...
      '''Quintic B''      C4\n',...
      '''Gaussian''       C infty\n'])
end
