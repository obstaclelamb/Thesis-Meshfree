function plotCompleteness(str,RK,varargin)
  x = RK.X;
  f = RK.f;
  C = RK.C;
  dC= RK.DC;
  ddC=RK.DDC;
  
  syms x_sym;
  f_sym = f(x_sym);
  df = diff(f_sym,x_sym);
  ddf = diff(df,x_sym);
  y_f = f(x);
  y_df = eval(subs(df,x_sym,x));
  y_ddf = eval(subs(ddf,x_sym,x));
  
  figure;  hold on;
  subplot(3,1,1)
  plot(x,C,'or',x,y_f,'k-')
  title([str,' completeness of N0'])
  legend('RK','exact')
  subplot(3,1,2)
  plot(x,dC,'or',x,y_df,'k-')
  title([str,' completeness of N1'])
  legend('RK','exact')
  subplot(3,1,3)
  plot(x,ddC,'or',x,y_ddf,'k-')
  title([str,' completeness of N2'])
  legend('RK','exact')

  str = ['DEG: ',num2str(RK.DEG),', kernel: ',RK.kernel]; 
  sgtitle(str)
end