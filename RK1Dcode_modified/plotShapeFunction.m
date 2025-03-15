function plotShapeFunction(str,N,x)
  figure;  hold on
  for k=1:length(N(1,:))
    Ntemp = N(:,k);
    plot(x,Ntemp,'LineWidth',2)
  end
  xlabel('x');
  title(str,'Interpreter','latex')
end